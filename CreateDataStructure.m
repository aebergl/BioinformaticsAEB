function DATA = CreateDataStructure(nRow,nCol,numIdColumns,numHeaderRows)
% DATA = CreateDataStructure(nRow,nCol,numIdColumns,numHeaderRows)
%
%   Initate DATA structure
%


DATA.Title = '';
DATA.ShortTitle = '';
DATA.X = zeros(nRow,nCol);
DATA.nRow = nRow;
DATA.nCol = nCol;
DATA.ColId = cell(DATA.nCol,1);
DATA.RowId = cell(DATA.nRow,1);
if numIdColumns > 1
    DATA.RowAnnotationFields = cell(numIdColumns-1,1);
    DATA.RowAnnotation = cell(DATA.nRow,numIdColumns-1);
else
    DATA.RowAnnotationFields = [];
    DATA.RowAnnotation = [];
end
if numHeaderRows > 1
    DATA.ColAnnotationFields = cell(numHeaderRows-1,1);
    DATA.ColAnnotation = cell(DATA.nCol,numHeaderRows-1);
else
    DATA.ColAnnotationFields = [];
    DATA.ColAnnotation = [];
end
DATA.Info.Source = '';
DATA.Info.Type = '';
DATA.Info.Platform = '';
DATA.Info.CreateDate = date;


end