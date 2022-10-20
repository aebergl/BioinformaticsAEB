function DATA = RemoveMaskCpGprobes(DATA,ChipType)


% Read masked information based on hg38
switch lower(ChipType)
    case {'450','m450','450k','m450k'}
        MaskData=readcell('/Users/bergluae/AEBERGL/DATA/METHYLATION/MASKED/HM450_hg38_Masked.txt','NumHeaderLines',1);
    case {'epic','850','850k'}
        MaskData=readcell('/Users/bergluae/AEBERGL/DATA/METHYLATION/MASKED/EPIC_hg38_Masked.txt','NumHeaderLines',1);
end

probes_MASKED = MaskData(strcmp('TRUE',MaskData(:,2)),1);
[~,~,DATA_indx] = intersect(probes_MASKED,DATA.ColId,'Stable');

DATA.X(:,DATA_indx) = [];

DATA.ColId(DATA_indx) = [];
DATA.nCol = size(DATA.X,2);
if ~isempty(DATA.ColAnnotation)
    DATA.ColAnnotation(DATA_indx,:) = [];
end

