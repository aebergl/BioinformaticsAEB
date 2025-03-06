function DATA  = AddIlluminaManifestAnnotation(DATA,AnnotationFile,varargin)



L = readlines('AnnotationFile');
AssayPos    = find(contains(a,'[Assay]'));
ControlsPos = find(contains(a,'[Controls]'));

ColumnIds = strsplit(L(AssayPos+1),',');



% Select black with annotation
L = L(AssayPos+2:ControlsPos-1);