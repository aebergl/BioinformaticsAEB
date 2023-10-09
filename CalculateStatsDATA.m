function RESULTS_DATA = CalculateStatsDATA(DATA,Type)


VarNames = ["n Values" "n MV" "Mean" "Median" "std" "Min" "Max" "Range" "IQR" "Skewness" "Kurtosis" "BC"]';
nVar = length(VarNames);

switch lower(Type)
    case 'samples'
        X = DATA.X';
        nRow = DATA.nRow;
        Title_txt = "Sample Statistics";
        RowId = DATA.RowId;
    case 'variables'
        X = DATA.X;
        nRow = DATA.nCol;
        Title_txt = "Variable Statistics";
        RowId = DATA.ColId;
    otherwise
        error('Type must either be ''samples'' or ''variables'' ')

end

nVal_X     = sum(~isnan(X),1);
nMV_X      = sum(isnan(X),1);
mean_X     = mean(X,1,'omitnan');
median_X   = median(X,1,'omitnan');
std_X      = std(X,1,'omitnan');
min_X      = min(X,1,'omitnan');
max_X      = max(X,1,'omitnan');
range_X    = range(X,1);
iqr_X      = iqr(X,1);
skewness_x = skewness(X,0,1 );
kurtosis_x = kurtosis(X,0,1) - 3;
BC_X       = ( (skewness_x.^2) + 1 ) ./ ( kurtosis_x + (3 .* (nVal_X-1) .^ 2 ./ ((nVal_X - 2) .* (nVal_X - 3))) );



RESULTS_DATA = CreateDataStructure(nRow,nVar,[],[]);


% Add Info
RESULTS_DATA.Title = Title_txt;
RESULTS_DATA.Info.Source = inputname(1);

RESULTS_DATA.ColId=VarNames;
RESULTS_DATA.RowId = DATA.RowId;
RESULTS_DATA.RowAnnotationFields = DATA.RowAnnotationFields;
RESULTS_DATA.RowAnnotation = DATA.RowAnnotation;



