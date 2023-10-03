function RESULTS_DATA = CalculateStatsDATA(DATA,Type)

nVar = 10;

switch lower(Type)
    case 'samples'
        X = DATA.X';
        nRow = DATA.nRow;
    case 'variables'
        X = DATA.X;
        nRow = DATA.nCol;
    otherwise
        error('Type must either be ''samples'' or ''variables'' ')

end

mean_X = mean(X,1,'omitnan');
median_X = median(X,1,'omitnan');
std_X  = std(X,1,'omitnan');
min_X  = min(X,1,'omitnan');
max_X  = max(X,1,'omitnan');
range_X = range(X,1);

RESULTS_DATA = CreateDataStructure(DATA.nRow,23,[],[]);