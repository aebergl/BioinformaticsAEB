function DATA = EditSamplesGroupDATA(DATA,RowAnnotationField,GroupsToRemove,KeepRemove)
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


if ~iscell(GroupsToRemove)
    GroupsToRemove = {GroupsToRemove};
end

indx_Row = strcmpi(RowAnnotationField,DATA.RowAnnotationFields);
if ~any(indx_Row)
    error('Could not find the given RowAnnotationField')
end

SampleIndx = false(DATA.nRow,1);

for i=1:numel(GroupsToRemove)
    SampleIndx = SampleIndx | strcmp(GroupsToRemove{i},DATA.RowAnnotation(:,indx_Row));
end

switch lower(KeepRemove)
    case 'keep'
        DATA  = EditSamplesDATA(DATA,DATA.RowId(SampleIndx),'Keep');        
    case 'remove'
        DATA  = EditSamplesDATA(DATA,DATA.RowId(SampleIndx),'Remove'); 
        
end

end




