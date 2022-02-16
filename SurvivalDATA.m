function [p_LR,fH,stats] = SurvivalDATA(DATA,SurvType,VarId,AnnoId,varargin)

GroupsToUse = [];

i=0;
while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'GroupsToUse')
        i = i + 1;
        GroupsToUse = varargin{i};
    end
end
if ~iscell(VarId)
    VarId = {VarId};
end
VarIndx = strcmp(VarId,DATA.ColId);


% which samples have survival data




end