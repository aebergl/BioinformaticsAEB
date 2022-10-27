function DATA = AddCpGProbeAnnotation(DATA,ProbeAnnotation)


ProbeId_Column = strcmp('IlmnID',ProbeAnnotation.ProbeAnnotationColumns);
if ~any(ProbeId_Column)
    ProbeId_Column = strcmp('TargetId',ProbeAnnotation.ProbeAnnotationColumns);
end


[~,Data_indx,Annotation_indx]=intersect(DATA.ColId,ProbeAnnotation.ProbeAnnotation(:,ProbeId_Column),'stable');

DATA.X = DATA.X(:,Data_indx);
DATA.nCol = size(DATA.X,2);
DATA.ColId = DATA.ColId(Data_indx);
DATA.ColAnnotation = ProbeAnnotation.ProbeAnnotation(Annotation_indx,:);
DATA.ColAnnotationFields = ProbeAnnotation.ProbeAnnotationColumns;


