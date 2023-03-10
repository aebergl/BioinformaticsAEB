function DATAm = Data2OldMet(DATA)

DATAm.X = DATA.X';
DATAm.P = zeros(size(DATAm.X));
DATAm.ProbeId = cellstr(DATA.ColId);
%DATAm.ProbeAnnotation = cellstr(DATA.ColAnnotation);
DATAm.ProbeAnnotation = (DATA.ColAnnotation);
% DATAm.ProbeAnnotationColumns = cellstr(DATA.ColAnnotationFields);
DATAm.SampleIds = cellstr(DATA.RowId);
DATAm.ProbeAnnotationColumns = (DATA.ColAnnotationFields);
%DATAm.SampleIds = (DATA.RowId);
DATAm.NumProbes = DATA.nCol;
DATAm.NumSamples = DATA.nRow;
% DATAm.SampleAnnotation = cellstr(DATA.RowAnnotation);
DATAm.SampleAnnotation = (DATA.RowAnnotation);
DATAm.SampleAnnotationColumns = DATA.RowAnnotationFields;