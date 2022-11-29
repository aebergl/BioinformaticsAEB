function RESULTS_DATA = CalculateDiffResults(DATA,Group1,Group2,GroupId,DataType,varargin)

MinNumSampleSize = 5;

if isempty(DataType)
    DataType = 'Beta-value';
end

if ~iscell(Group1)
    Group1 = {Group1};
end
if ~iscell(Group2)
    Group2 = {Group2};
end

GroupIdColumn = (strcmp(GroupId,DATA.RowAnnotationFields));
GroupIdColumn = GroupIdColumn & (cumsum(GroupIdColumn) == 1); % Picks the first occurence

SampleIndxMatrix = false(DATA.nRow,numel(Group1));
for i=1:numel(Group1)
    SampleIndxMatrix(:,i) = strcmp(Group1{i},DATA.RowAnnotation(:,GroupIdColumn));
end
SampleIndx_x = any(SampleIndxMatrix,2);

SampleIndxMatrix = false(DATA.nRow,numel(Group2));
for i=1:numel(Group2)
    SampleIndxMatrix(:,i) = strcmp(Group2{i},DATA.RowAnnotation(:,GroupIdColumn));
end
SampleIndx_y = any(SampleIndxMatrix,2);

SampleId_x = DATA.RowId(SampleIndx_x);
SampleId_y = DATA.RowId(SampleIndx_y);


RESULTS_DATA = CreateDataStructure(DATA.nCol,10,[],[]);

% Add Info
RESULTS_DATA.Title = 'Difference analysis results';
RESULTS_DATA.Info.Source = inputname(1);
RESULTS_DATA.Info.Platform = DataType;
RESULTS_DATA.Info.GroupID = GroupId;
RESULTS_DATA.Info.Group1 = Group1;
RESULTS_DATA.Info.Group2 = Group2;
RESULTS_DATA.Info.numSamples_x = length(SampleId_x);
RESULTS_DATA.Info.numSamples_y = length(SampleId_y);
RESULTS_DATA.Info.SampleId_x = SampleId_x;
RESULTS_DATA.Info.SampleId_y = SampleId_y;
RESULTS_DATA.Info.MinNumSampleSize = MinNumSampleSize;
RESULTS_DATA.RowId = DATA.ColId;
RESULTS_DATA.RowAnnotation = DATA.ColAnnotation;
RESULTS_DATA.RowAnnotationFields = DATA.ColAnnotationFields;



VarNames = {'p t-test','q t-test','fdr t-test','Delta Average','Average x','Average y','p Bartlett','q Bartlett','fdr Bartlett','Min Samples'}';
RESULTS_DATA.ColId=VarNames;

RESULTS_DATA.ColAnnotationFields = {'VarId'};
RESULTS_DATA.ColAnnotation = RESULTS_DATA.ColId;



p_TT = ones(DATA.nCol,1) * NaN;
DeltaAverage = zeros(DATA.nCol,1) * NaN;
Average_x = zeros(DATA.nCol,1) * NaN;
Average_y = zeros(DATA.nCol,1) * NaN;
p_Bartlett = ones(DATA.nCol,1) * NaN;
MinNum = zeros(DATA.nCol,1) * NaN;

%for i=1:1000%DATA.NumProbes
X_x = DATA.X(SampleIndx_x,:);
X_y = DATA.X(SampleIndx_y,:);
if strcmpi('M-value',DataType)
    X_x = B2M(X_x);
    X_y = B2M(X_y);
end


parfor i=1:DATA.nCol
    x = X_x(:,i);
    y = X_y(:,i);
    x(isnan(x)) = [];
    y(isnan(y)) = [];
    MinNum(i) = min([length(x),length(y)]);
    if MinNum(i) > MinNumSampleSize
        [~,p_TT(i)] = ttest2(x,y,0.05,'both','unequal');
        Average_x(i) = mean(x,'omitnan');
        Average_y(i) = mean(y,'omitnan');
        p_Bartlett(i)=vartestn([x;y],[ones(size(x));ones(size(y))*2],'Display','off','TestType','Bartlett');
        DeltaAverage(i) = Average_y(i) - Average_x(i);
    end
end


RESULTS_DATA.X(:,1) = p_TT;
RESULTS_DATA.X(:,4) = DeltaAverage;
RESULTS_DATA.X(:,5) = Average_x;
RESULTS_DATA.X(:,6) = Average_y;
RESULTS_DATA.X(:,7) = p_Bartlett;
RESULTS_DATA.X(:,10) = MinNum;

try
    [~, RESULTS_DATA.X(:,2),~] = mafdr(p_TT);
catch
    RESULTS_DATA.X(:,2) = ones(DATA.nCol,1);
end
try
    [RESULTS_DATA.X(:,3)] = mafdr(p_TT,'BHFDR',true);
catch
    RESULTS_DATA.X(:,3) = ones(DATA.nCol,1);
end

try
    [~, RESULTS_DATA.X(:,8),~] = mafdr(p_Bartlett);
catch
    RESULTS_DATA.X(:,8) = ones(DATA.nCol,1);
end
try
    [RESULTS_DATA.X(:,9)] = mafdr(p_Bartlett,'BHFDR',true);
catch
    RESULTS_DATA.X(:,9) = ones(DATA.nCol,1);
end

