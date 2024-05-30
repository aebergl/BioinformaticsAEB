function DATA = ConvertDataOld2New(DATA)


DATA = renameStructField(DATA, 'NumSamples', 'nRow');
DATA = renameStructField(DATA, 'NumProbes', 'nCol');

DATA = renameStructField(DATA, 'ProbeId', 'ColId');
DATA = renameStructField(DATA, 'SampleId', 'RowId');

DATA = renameStructField(DATA, 'SampleAnnotationFields', 'RowAnnotationFields');
DATA = renameStructField(DATA, 'ProbeAnnotationFields', 'ColAnnotationFields');

DATA = renameStructField(DATA, 'SampleAnnotation', 'RowAnnotation');
DATA = renameStructField(DATA, 'ProbeAnnotation', 'ColAnnotation');

DATA = rmfield(DATA,'ProbeXDim');
DATA = rmfield(DATA,'SampleXDim');

