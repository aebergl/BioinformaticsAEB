function DATA = CreateAEBDataStructure(nRow,nCol,numIdColumns,numHeaderRows)

DATA.Title = '';
DATA.ShortTitle = '';
DATA.X = zeros(nRow,nCol);
DATA.nRow = nRow;
DATA.nCol = nCol;
DATA.ColId = cell(DATA.nCol,1);
DATA.RowId = cell(DATA.nRow,1);
DATA.RowAnnotationFields = cell(numIdColumns,1);
DATA.RowAnnotation = cell(DATA.nRow,numIdColumns);
DATA.ColAnnotationFields = cell(numHeaderRows,1);
DATA.ColAnnotation = cell(DATA.nCol,numHeaderRows);
DATA.Info.Source = '';
DATA.Info.Type = '';
DATA.Info.Platform = '';
DATA.Info.CreateDate = date;


end