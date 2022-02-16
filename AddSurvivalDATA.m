function DATA = AddSurvivalDATA(DATA,IdToUse,SURVIVAL)
% USAGE:
%   DATA = AddSurvivalDATA(DATA,SURVIVAL) Add Survival structure to DATA
%   structure
%
% INPUTS:
% * DATA: DATA structure
% * IdToUse: Row identifier to use, [], RowId is used
% * SURVIVAL: SURVIVAL structure 
%
% OUTPUTS:
% * DATA: DATA structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Anders Berglund, 2020 aebergl at gmail.com                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isempty(IdToUse) || strcmpi(IdToUse,'RowId')
    DATAid = DATA.RowId;
else
    indx = strcmpi(IdToUse,DATA.RowAnnotationFields);
    if any(indx)
        DATAid = DATA.RowAnnotation(:,indx);
    else
        error('Could not find the given row identifier')
    end   
end

nCharDATAid = length(DATAid{1});

SURVid = SURVIVAL.RowId;
nCharSURVid = length(SURVid{1});
nSURV_Items = length(SURVIVAL.SurvivalTypes);

% Check if the IDs have the same length
if nCharDATAid < nCharSURVid
    SURVid = cellfun(@(x) x(1:nCharDATAid), SURVid, 'UniformOutput', false);
elseif nCharDATAid > nCharSURVid
    DATAid = cellfun(@(x) x(1:nCharSURVid), DATAid, 'UniformOutput', false);
end

%Create DATA SURVIVAL Object
DATA.SURVIVAL.RowId = cell(DATA.nRow,1);
DATA.SURVIVAL.RowId(:) = {''};
DATA.SURVIVAL.SurvivalTypes = SURVIVAL.SurvivalTypes;
DATA.SURVIVAL.Units = SURVIVAL.Units;
DATA.SURVIVAL.SurvEvent = cell(DATA.nRow,nSURV_Items);
DATA.SURVIVAL.SurvEvent(:) = {'NA'};
DATA.SURVIVAL.SurvTime = ones(DATA.nRow,nSURV_Items) * NaN;

% Find matching Ids
[indx_DATA, location_SURV]  = ismember(DATAid, SURVid);
location_SURV = location_SURV(location_SURV>0);
DATA.SURVIVAL.RowId(indx_DATA) = SURVIVAL.RowId(location_SURV);
DATA.SURVIVAL.SurvEvent(indx_DATA,:) = SURVIVAL.SurvEvent(location_SURV,:);
DATA.SURVIVAL.SurvTime(indx_DATA,:) = SURVIVAL.SurvTime(location_SURV,:);






    