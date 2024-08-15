function DATA_OUT = CalculateDiffVal_DATA(DATA, RowIds1, RowIds2, varargin)

SampleIdentifier = false;

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


indx1 = ismember(DATA_ID,RowIds1);
indx2 = ismember(DATA_ID,RowIds2);

if sum(indx1) ~= sum(indx2)
    error('Different number of rows found')
end

if any(all([indx1 indx2],2))
    error('Same row found in both RowIds1 and RowIds2')
end

X1 = DATA.X(indx1,:);
X2 = DATA.X(indx2,:);

Xdiff = X2 - X1;

DATA_OUT = CreateDataStructure(sum(indx1),DATA.nCol,[],[]);
DATA_OUT.X = Xdiff;
DATA_OUT.RowAnnotationFields = DATA.RowAnnotationFields;
DATA_OUT.RowAnnotation = DATA.RowAnnotation(indx1,:);
DATA_OUT.ColId = strcat(DATA.ColId," change");
DATA_OUT.RowId = strcat(DATA_ID(indx2)," - ",DATA_ID(indx2));

if isfield(DATA,'SURVIVAL')
    DATA_OUT.SURVIVAL.RowId =  DATA.SURVIVAL.RowId(indx1);
    DATA_OUT.SURVIVAL.SurvivalTypes =  DATA.SURVIVAL.SurvivalTypes;
    DATA_OUT.SURVIVAL.Units =  DATA.SURVIVAL.Units;
    DATA_OUT.SURVIVAL.SurvEvent =  DATA.SURVIVAL.SurvEvent(indx1,:);
    DATA_OUT.SURVIVAL.SurvTime =  DATA.SURVIVAL.SurvTime(indx1,:);
end