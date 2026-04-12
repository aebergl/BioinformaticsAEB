function DATA_v1v2 = MergeEPICv1v2(DATA_v1,DATA_v2)

load("/Users/berglund.anders/Documents/DATA/ANNOTATION/METHYLATION/EPICv2/EPICv2EPICv1Mapping/EPICv1EPICv2_Mapping_strict.mat",'EPICv1EPICv2_Mapping_strict')

% Remove EPICv1 missing CPGs
Lia = ismember(EPICv1EPICv2_Mapping_strict(:,2),DATA_v1.ColId);
EPICv1EPICv2_Mapping_strict = EPICv1EPICv2_Mapping_strict(Lia,:);

% Remove EPICv2 missing CPGs
Lia = ismember(EPICv1EPICv2_Mapping_strict(:,1),DATA_v2.ColId);
EPICv1EPICv2_Mapping_strict = EPICv1EPICv2_Mapping_strict(Lia,:);

[~,~,indx_v1] = intersect(EPICv1EPICv2_Mapping_strict(:,2),DATA_v1.ColId,'stable');
[~,~,indx_v2] = intersect(EPICv1EPICv2_Mapping_strict(:,1),DATA_v2.ColId,'stable');

DATA_v1v2 = CreateDataStructure(DATA_v1.nRow+DATA_v2.nRow,length(indx_v1),[],[]);
DATA_v1v2.RowId = [DATA_v1.RowId;DATA_v2.RowId];
DATA_v1v2.X = [DATA_v1.X(:,indx_v1);DATA_v2.X(:,indx_v2)];
DATA_v1v2.ColId = DATA_v2.ColId(indx_v2);
DATA_v1v2.ColAnnotation = DATA_v2.ColAnnotation(indx_v2,:);
DATA_v1v2.ColAnnotationFields = DATA_v2.ColAnnotationFields;
DATA_v1v2.RowAnnotation = [repmat("EPICv1",DATA_v1.nRow,1);repmat("EPICv2",DATA_v2.nRow,1)];
DATA_v1v2.RowAnnotationFields = "ChipType";
