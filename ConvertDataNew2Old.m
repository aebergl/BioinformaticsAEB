function DATA = ConvertDataNew2Old(DATA)


DATA = renameStructField(DATA, 'nRow', 'NumSamples');
DATA = renameStructField(DATA, 'nCol', 'NumProbes');

DATA = renameStructField(DATA, 'ColId', 'ProbeId');
DATA = renameStructField(DATA, 'RowId', 'SampleId');

DATA = renameStructField(DATA, 'RowAnnotationFields', 'SampleAnnotationFields');
DATA = renameStructField(DATA, 'ColAnnotationFields', 'ProbeAnnotationFields');

DATA = renameStructField(DATA, 'RowAnnotation', 'SampleAnnotation');
DATA = renameStructField(DATA, 'ColAnnotation', 'ProbeAnnotation');


DATA.ProbeXDim = 2;
DATA.SampleXDim = 1;
