function RESULTS_DATA = CalculateCorrelationResults(DATA,Y,Group1,GroupId,DataType,varargin)


MinNumSampleSize = 5;

if isempty(DataType)
    DataType = 'Beta-value';
end

if ~iscell(Group1)
    Group1 = {Group1};
end


if ~isempty(GroupId)
    GroupIdColumn = (strcmp(GroupId,DATA.RowAnnotationFields));
    GroupIdColumn = GroupIdColumn & (cumsum(GroupIdColumn) == 1); % Picks the first occurence

    SampleIndxMatrix = false(DATA.nRow,numel(Group1));
    for i=1:numel(Group1)
        SampleIndxMatrix(:,i) = strcmp(Group1{i},DATA.RowAnnotation(:,GroupIdColumn));
    end
    SampleIndx_x = any(SampleIndxMatrix,2);
    SampleId_x = DATA.RowId(SampleIndx_x);
else
    SampleIndx_x = true(DATA.nRow,1);
    SampleId_x = DATA.RowId;
    GroupId = 'All samples';
    Group1 = 'All samples';
end

RESULTS_DATA = CreateDataStructure(DATA.nCol,10,[],[]);

% Add Info
RESULTS_DATA.Title = 'Correlation analysis results';
RESULTS_DATA.Info.Source = inputname(1);

RESULTS_DATA.Info.DataType = DataType;
RESULTS_DATA.Info.GroupID = GroupId;
RESULTS_DATA.Info.Group1 = Group1;
RESULTS_DATA.Info.MinNumSampleSize = MinNumSampleSize;
RESULTS_DATA.Info.SampleId_x = SampleId_x;
RESULTS_DATA.Info.numSamples_x = length(SampleId_x);
RESULTS_DATA.RowId = DATA.ColId;
RESULTS_DATA.RowAnnotation = DATA.ColAnnotation;
RESULTS_DATA.RowAnnotationFields = DATA.ColAnnotationFields;

VarNames = {'r Pearson','p Pearson','q Pearson','fdr Pearson','r Spearman','p Spearman','q Spearman','fdr Spearman','Range','Num Values'}';
RESULTS_DATA.ColId=VarNames;

RESULTS_DATA.ColAnnotationFields = {'VarId'};
RESULTS_DATA.ColAnnotation = RESULTS_DATA.ColId;

X_x = DATA.X(SampleIndx_x,:);

Y = Y(SampleIndx_x);
if strcmpi('M-value',DataType)
    X_x = B2M(X_x);
end


r_Pearson = ones(DATA.nCol,1) * NaN;
p_Pearson = ones(DATA.nCol,1) * NaN;
r_Spearman = ones(DATA.nCol,1) * NaN;
p_Spearman = ones(DATA.nCol,1) * NaN;
RangeVal  = ones(DATA.nCol,1) * NaN;
nVal = ones(DATA.nCol,1) * NaN;


parfor i=1:DATA.nCol
    x = X_x(:,i);
    nVal(i) = sum(~isnan(x));
    if nVal(i) > MinNumSampleSize
        RangeVal(i) = range(x);
        [r_Pearson(i), p_Pearson(i)] = corr(x,Y,'Type','Pearson','Rows','Pairwise');
        try
            [r_Spearman(i), p_Spearman(i)] = corr(x,Y,'Type','Spearman','Rows','Pairwise');
        catch
            [r_Spearman(i)] = corr(x,Y,'Type','Spearman','Rows','Pairwise');
            p_Spearman(i) = 1e-10;
        end
    end
end
RESULTS_DATA.X(:,1) = r_Pearson;
RESULTS_DATA.X(:,2) = p_Pearson;
RESULTS_DATA.X(:,5) = r_Spearman;
RESULTS_DATA.X(:,6) = p_Spearman;
RESULTS_DATA.X(:,9) = RangeVal;
RESULTS_DATA.X(:,10) = nVal;
try
    [~, RESULTS_DATA.X(:,3),~] = mafdr(p_Pearson);
catch
    RESULTS_DATA.X(:,3) = ones(DATA.nCol,1);
end

try
    [RESULTS_DATA.X(:,4)] = mafdr(p_Pearson,'BHFDR',true);
catch
    RESULTS_DATA.X(:,4) = ones(DATA.nCol,1);
end


try
    [~, RESULTS_DATA.X(:,7),~] = mafdr(p_Spearman);
catch
    RESULTS_DATA.X(:,7) = ones(DATA.nCol,1);
end

try
    [RESULTS_DATA.X(:,8)] = mafdr(p_Spearman,'BHFDR',true);
catch
    RESULTS_DATA.X(:,8) = ones(DATA.nCol,1);
end
