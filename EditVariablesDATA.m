function DATA  = EditVariablesDATA(DATA,InputIds,KeepRemove,varargin)
%DATA  = EditSamplesDATA(DATA,InputIds,KeepRemove,varargin)
%
%   Remove/Keeps samples from DATA
%
% INPUT
%   DATA                DATA structure
%   InputIds            List of variables to keep or remove
%   KeepRemove          'Keep' or 'Remove'
%
%   options ---------------------------------------
%
%   'VariableIdentifier'    Sampleidentifier to be used, default DATA.SampleId
%   'Truncate'              Number of characters to use for matching, defualt all
% Anders Berglund


Truncate = 0;
VariableIdentifier = false;

i=0;
while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'VariableIdentifier')
        i = i + 1;
        VariableIdentifier = varargin{i};
    elseif strcmpi(varargin{i},'Truncate')
        i = i + 1;
        Truncate = varargin{i};
    end
end

if VariableIdentifier
    indx = strcmpi(VariableIdentifier,DATA.ColAnnotationFields);
    if ~any(indx)
        error('Error. \n%s not found in DATA.ColAnnotationFields',VariableIdentifier);
    elseif sum(indx) > 1
        error('Warning. \nMultiple matches for %s found in DATA.ColAnnotationFields',VariableIdentifier);
    else
        DATA_ID = DATA.ColAnnotation(:,indx);
    end
else
    DATA_ID = DATA.ColId;
end

if Truncate
    DATA_ID = cellfun(@(x) x(1:Truncate), DATA_ID, 'UniformOutput', false);
    InputIds = cellfun(@(x) x(1:Truncate), InputIds, 'UniformOutput', false);
end

IdUnique = unique(InputIds);
if length(IdUnique) < length(InputIds)
    warning('WARNING! Not all input IDs are unique in EditSamplesDATA')
end


[indx] = ismember(DATA_ID,InputIds);

if any(indx)
    switch lower(KeepRemove)
        case 'keep'
            DATA.X = DATA.X(:,indx);
            DATA.ColId = DATA.ColId(indx);
            DATA.nCol = size(DATA.X,2);
            if ~isempty(DATA.ColAnnotation)
                DATA.ColAnnotation = DATA.ColAnnotation(indx,:);
            end
        case 'remove'
            %Remove samples
            DATA.X(:,indx) = [];
            DATA.ColId(indx) = [];
            DATA.nCol = size(DATA.X,2);
            if ~isempty(DATA.ColAnnotation)
                DATA.ColAnnotation(indx,:) = [];
            end
    end

else
    warning('WARNING! No matching Ids found, returning original DATA')
end


end