function DATA  = EditVariablesDATA(DATA,InputIds,KeepRemove,varargin)
%DATA  = EditVariablesDATA(DATA,InputIds,KeepRemove,varargin)
%
%   Remove/Keeps variables from DATA
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
%   'Stable'                Returns the data in the order of InputIds
% Anders Berglund


Truncate = 0;
VariableIdentifier = false;
Stable = false;

i=0;
while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'VariableIdentifier')
        i = i + 1;
        VariableIdentifier = varargin{i};
    elseif strcmpi(varargin{i},'Truncate')
        i = i + 1;
        Truncate = varargin{i};
    elseif strcmpi(varargin{i},'Stable')
        Stable = true;
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


if Stable
    if length(IdUnique) == DATA.nCol && strcmpi('keep',KeepRemove) % Simple if there is unique list of IDs
        [~,~,indx] = intersect(InputIds,DATA_ID,'Stable');
        DATA.X = DATA.X(:,indx);
        DATA.RowId = DATA.ColId(indx);
        DATA.nRow = size(DATA.X,2);
        if ~isempty(DATA.ColAnnotation)
            DATA.ColAnnotation = DATA.ColAnnotation(indx,:);
        end
    else

        SampleIndx = false(DATA.nCol,length(IdUnique));
        for i=1:length(IdUnique)
            SampleIndx(:,i) = ismember(DATA_ID,InputIds(i)); % Find samples for each ID
        end
        sum_indx = sum(SampleIndx,2);
        if max(sum_indx) > 1
            error('Multiple samples found for differenb IDs')
        end
        [indx,~]=find(SampleIndx);
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


else
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

end