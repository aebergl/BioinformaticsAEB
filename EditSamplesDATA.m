function DATA  = EditSamplesDATA(DATA,InputIds,KeepRemove,varargin)
%DATA  = EditSamplesDATA(DATA,InputIds,KeepRemove,varargin)
%
%   Remove/Keeps samples from DATA
%
% INPUT
%   DATA                DATA structure
%   InputIds            List of samples to remove
%   KeepRemove          'Keep' or 'Remove'
%
%   options ---------------------------------------
%
%   'SampleIdentifier'      Sampleidentifier to be used, default DATA.SampleId
%   'Truncate'              Number of characters to use for matching, defualt all
% Anders Berglund


Truncate = false;
SampleIdentifier = false;
Stable = false;

i=0;
while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'SampleIdentifier')
        i = i + 1;
        SampleIdentifier = varargin{i};
    elseif strcmpi(varargin{i},'Truncate')
        i = i + 1;
        Truncate = varargin{i};
    elseif strcmpi(varargin{i},'Stable')
        Stable = true;
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

if Truncate
    DATA_ID = cellfun(@(x) x(1:Truncate), DATA_ID, 'UniformOutput', false);
    InputIds = cellfun(@(x) x(1:Truncate), InputIds, 'UniformOutput', false);
end

IdUnique = unique(InputIds,'stable');
if length(IdUnique) < length(InputIds)
    warning('WARNING! Not all input IDs are unique in EditSamplesDATA')
end

if Stable

    if length(IdUnique) == DATA.nRow && strcmpi('keep',KeepRemove)
        [~,~,indx] = intersect(InputIds,DATA_ID,'Stable');
        DATA.X = DATA.X(indx,:);
        DATA.RowId = DATA.RowId(indx);
        DATA.nRow = size(DATA.X,1);
        DATA.RowAnnotation = DATA.RowAnnotation(indx,:);
        if isfield(DATA,'SURVIVAL')
            DATA.SURVIVAL.RowId =  DATA.SURVIVAL.RowId(indx);
            DATA.SURVIVAL.SurvEvent =  DATA.SURVIVAL.SurvEvent(indx,:);
            DATA.SURVIVAL.SurvTime =  DATA.SURVIVAL.SurvTime(indx,:);
        end
    else

        SampleIndx = false(DATA.nRow,length(IdUnique));
        for i=1:length(IdUnique)
            SampleIndx(:,i) = ismember(DATA_ID,InputIds(i));
        end
        sum_indx = sum(SampleIndx,2);
        if max(sum_indx) > 1
            error('Multiple match samples found for differenb IDs')
        end
        [indx,~]=find(SampleIndx);
        if any(indx)
            switch lower(KeepRemove)
                case 'keep'
                    DATA.X = DATA.X(indx,:);
                    DATA.RowId = DATA.RowId(indx);
                    DATA.nRow = size(DATA.X,1);
                    DATA.RowAnnotation = DATA.RowAnnotation(indx,:);
                    if isfield(DATA,'SURVIVAL')
                        DATA.SURVIVAL.RowId =  DATA.SURVIVAL.RowId(indx);
                        DATA.SURVIVAL.SurvEvent =  DATA.SURVIVAL.SurvEvent(indx,:);
                        DATA.SURVIVAL.SurvTime =  DATA.SURVIVAL.SurvTime(indx,:);
                    end
                case 'remove'
                    %Remove samples
                    DATA.X(indx,:) = [];
                    DATA.RowId(indx) = [];
                    DATA.nRow = size(DATA.X,1);
                    DATA.RowAnnotation(indx,:) = [];
                    if isfield(DATA,'SURVIVAL')
                        DATA.SURVIVAL.RowId(indx) = [];
                        DATA.SURVIVAL.SurvEvent(indx,:) =  [];
                        DATA.SURVIVAL.SurvTime(indx,:) =  [];
                    end
            end

        else
            warning('WARNING! No matching Ids found, returning original DATA')
        end
    end

else
    indx = ismember(DATA_ID,InputIds);

    if any(indx)
        switch lower(KeepRemove)
            case 'keep'
                DATA.X = DATA.X(indx,:);
                DATA.RowId = DATA.RowId(indx);
                DATA.nRow = size(DATA.X,1);
                if ~isempty(DATA.RowAnnotation )
                    DATA.RowAnnotation = DATA.RowAnnotation(indx,:);
                end
                if isfield(DATA,'SURVIVAL')
                    DATA.SURVIVAL.RowId =  DATA.SURVIVAL.RowId(indx);
                    DATA.SURVIVAL.SurvEvent =  DATA.SURVIVAL.SurvEvent(indx,:);
                    DATA.SURVIVAL.SurvTime =  DATA.SURVIVAL.SurvTime(indx,:);
                end
            case 'remove'
                %Remove samples
                DATA.X(indx,:) = [];
                DATA.RowId(indx) = [];
                DATA.nRow = size(DATA.X,1);
                if ~isempty(DATA.RowAnnotation )
                    DATA.RowAnnotation(indx,:) = [];
                end
                if isfield(DATA,'SURVIVAL')
                    DATA.SURVIVAL.RowId(indx) = [];
                    DATA.SURVIVAL.SurvEvent(indx,:) =  [];
                    DATA.SURVIVAL.SurvTime(indx,:) =  [];
                end
        end

    else
        warning('WARNING! No matching Ids found, returning original DATA')
    end
end


end