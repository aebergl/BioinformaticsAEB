function DATA = TransposeDataAEB(DATA)

DATA.X      = DATA.X';
DATA.nRow   = size(DATA.X,1);
DATA.nCol   = size(DATA.X,2);

tmp         = DATA.ColId;
DATA.ColId  = DATA.RowId;
DATA.RowId  = tmp;

tmp                         = DATA.RowAnnotationFields;
DATA.RowAnnotationFields    = DATA.ColAnnotationFields;
DATA.ColAnnotationFields    = tmp;

tmp                         = DATA.RowAnnotation;
DATA.RowAnnotation          = DATA.ColAnnotation;
DATA.ColAnnotation          = tmp;


