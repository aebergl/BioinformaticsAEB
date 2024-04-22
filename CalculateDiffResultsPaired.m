function RESULTS_DATA = CalculateDiffResultsPaired(DATA,SampleId_x,SampleId_y,GroupId,DataType,varargin)

MinNumSampleSize = 2;

if isempty(DataType)
    DataType = 'Beta-value';
end

GroupIdColumn = (strcmp(GroupId,DATA.RowAnnotationFields));
GroupIdColumn = GroupIdColumn & (cumsum(GroupIdColumn) == 1); % Picks the first occurence

[~,SampleIndx_x,~] = intersect(DATA.RowAnnotation(:,GroupIdColumn),SampleId_x,'stable');
[~,SampleIndx_y,~] = intersect(DATA.RowAnnotation(:,GroupIdColumn),SampleId_y,'stable');

if length(SampleIndx_x) ~= length(SampleIndx_y)
    error('There need to my same number of samples in each group')
end


RESULTS_DATA = CreateDataStructure(DATA.nCol,7,[],[]);

% Add Info
RESULTS_DATA.Title = 'Paired Difference analysis results';
RESULTS_DATA.Info.Source = inputname(1);
RESULTS_DATA.Info.DataType = DataType;
RESULTS_DATA.Info.GroupID = GroupId;
RESULTS_DATA.Info.numSamples_x = length(SampleId_x);
RESULTS_DATA.Info.numSamples_y = length(SampleId_y);
RESULTS_DATA.Info.SampleId_x = SampleId_x;
RESULTS_DATA.Info.SampleId_y = SampleId_y;
RESULTS_DATA.Info.MinNumSampleSize = MinNumSampleSize;
RESULTS_DATA.RowId = DATA.ColId;
RESULTS_DATA.RowAnnotation = DATA.ColAnnotation;
RESULTS_DATA.RowAnnotationFields = DATA.ColAnnotationFields;



VarNames = {'p t-test','q t-test','fdr t-test','Delta Average','Average A','Average B','Min Samples'}';
RESULTS_DATA.ColId=VarNames;

RESULTS_DATA.ColAnnotationFields = {'VarId'};
RESULTS_DATA.ColAnnotation = RESULTS_DATA.ColId;


p_TT = ones(DATA.nCol,1) * NaN;
DeltaAverage = zeros(DATA.nCol,1) * NaN;
Average_x = zeros(DATA.nCol,1) * NaN;
Average_y = zeros(DATA.nCol,1) * NaN;
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
    indx_mv = isnan(x) | isnan(y);
    x(indx_mv) = [];
    y(indx_mv) = [];
    MinNum(i) = min([length(x),length(y)]);
    if MinNum(i) > MinNumSampleSize
        [~,p_TT(i)] = ttest2(x,y,0.05,'both');
        Average_x(i) = mean(x,'omitnan');
        Average_y(i) = mean(y,'omitnan');
        %DeltaAverage(i) = Average_y(i) - Average_x(i);
        DeltaAverage(i) = mean(y - x);
    end
end

RESULTS_DATA.X(:,1) = p_TT;
RESULTS_DATA.X(:,4) = DeltaAverage;
RESULTS_DATA.X(:,5) = Average_x;
RESULTS_DATA.X(:,6) = Average_y;
RESULTS_DATA.X(:,7) = MinNum;

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


