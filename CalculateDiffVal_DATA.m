function DATA_OUT = CalculateDiffVal_DATA(DATA, RowIdsNormal, RowIdsTumor, varargin)

SampleIdentifier = false;
KeepTumorId = true;
i=0;
while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'SampleIdentifier')
        i = i + 1;
        SampleIdentifier = varargin{i};

    end
end

if SampleIdentifier
    indx = strcmpi(SampleIdentifier,DATA.RowAnnotationFields);
    if ~any(indx)
        error('Error. \n%s not found in DATA.RowAnnotationFields',SampleIdentifier);
    elseif sum(indx) > 1
        error('Warning. \nMultiple matches for %s found in DATA.RowAnnotationFields',SampleIdentifier);
    else
        DATA_ID = DATA.RowAnnotation(:,indx);
    end
else
    DATA_ID = DATA.RowId;
end


[~,indxNormal] = ismember(RowIdsNormal,DATA_ID);
[~,indxTumor] = ismember(RowIdsTumor,DATA_ID);

if length(indxNormal) ~= length(indxTumor)
    error('Different number of rows found')
end

SameId = intersect(indxNormal, indxTumor);
if ~isempty(SameId)
    error('Same row found in both RowIds1 and RowIds2')
end

XNormal = DATA.X(indxNormal,:);
XTumor = DATA.X(indxTumor,:);

Xdiff = XTumor - XNormal;

DATA_OUT = CreateDataStructure(length(indxNormal),DATA.nCol,[],[]);
DATA_OUT.X = Xdiff;
DATA_OUT.RowAnnotationFields = DATA.RowAnnotationFields;
DATA_OUT.RowAnnotation = DATA.RowAnnotation(indxTumor,:);

if KeepTumorId
    DATA_OUT.ColId = DATA.ColId;
    DATA_OUT.RowId = DATA_ID(indxTumor);
    DATA_OUT.ColAnnotationFields = DATA.ColAnnotationFields;
    DATA_OUT.ColAnnotation = DATA.ColAnnotation;

else
    DATA_OUT.ColId = strcat(DATA.ColId," change");
    DATA_OUT.RowId = strcat(DATA_ID(indxTumor)," - ",DATA_ID(indxTumor));
end

if isfield(DATA,'SURVIVAL')
    DATA_OUT.SURVIVAL.RowId =  DATA.SURVIVAL.RowId(indxNormal);
    DATA_OUT.SURVIVAL.SurvivalTypes =  DATA.SURVIVAL.SurvivalTypes;
    DATA_OUT.SURVIVAL.Units =  DATA.SURVIVAL.Units;
    DATA_OUT.SURVIVAL.SurvEvent =  DATA.SURVIVAL.SurvEvent(indxNormal,:);
    DATA_OUT.SURVIVAL.SurvTime =  DATA.SURVIVAL.SurvTime(indxNormal,:);
end