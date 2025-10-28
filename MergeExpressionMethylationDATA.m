function [DATA_E, DATA_M] = MergeExpressionMethylationDATA(DATA_E,GeneIdE,SampleIdColumnE,DATA_M,SampleIdColumnM,varargin)

VariableIdentifier = [];

GeneId = ConvertStr(GeneIdE,'string');
SampleIdColumnE = ConvertStr(SampleIdColumnE,'string');
SampleIdColumnM = ConvertStr(SampleIdColumnM,'string');

Truncate = 16;


if nargin > 5
    ArgsList = {'VariableIdentifier'};
    for j=1:2:numel(varargin)

        ArgType = varargin{j};
        ArgVal = varargin{j+1};
        if ~strncmpi(ArgType,ArgsList,numel(ArgType))
            error('Invalid input option: %s', ArgType);
        else
            switch lower(ArgType)

                case 'variableidentifier'
                    VariableIdentifier = ArgVal;
                    VariableIdentifier = ConvertStr(VariableIdentifier,'string');
            end
        end
    end
end

% Selection of Y variable
if ~isempty(VariableIdentifier)
    indx_VarId = strcmpi(VariableIdentifier,DATA_E.ColAnnotationFields);
    if ~any(indx_VarId)
        error('Error. \n%s not found in DATA_E.ColAnnotationFields',VariableIdentifier);
    elseif sum(indx_VarId) > 1
        error('Warning. \nMultiple matches for %s found in DATA_E.ColAnnotationFields',VariableIdentifier);
    else
        DATA_ID = DATA_E.ColAnnotation(:,indx_VarId);
    end
else
    DATA_ID = DATA_E.ColId;
end
indx = ismember(DATA_ID,GeneIdE);
if sum(indx) > 1
    error('Multiple entries found for %s',GeneIdE);
else
    DATA_E=EditVariablesDATA(DATA_E,DATA_E.ColId(indx),'Keep');
end


if ~isempty(SampleIdColumnE)
    indx = strcmpi(SampleIdColumnE,DATA_E.RowAnnotationFields);
    if ~any(indx)
        error('Error. \n%s not found in DATA_E.RowAnnotationFields',SampleIdColumnE);
    elseif sum(indx) > 1
        error('Warning. \nMultiple matches for %s found in DATA_E.RowAnnotationFields',SampleIdColumnE);
    else
        DATA_ID_E = DATA_E.RowAnnotation(:,indx);
    end
else
    DATA_ID_E = DATA_E.RowId;
end

if ~isempty(SampleIdColumnE)
    indx = strcmpi(SampleIdColumnE,DATA_M.RowAnnotationFields);
    if ~any(indx)
        error('Error. \n%s not found in DATA_M.RowAnnotationFields',SampleIdColumnM);
    elseif sum(indx) > 1
        error('Warning. \nMultiple matches for %s found in DATA_M.RowAnnotationFields',SampleIdColumnM);
    else
        DATA_ID_M = DATA_M.RowAnnotation(:,indx);
    end
else
    DATA_ID_M = DATA_M.RowId;
end



if Truncate
    DATA_ID_E = cellfun(@(x) x(1:Truncate), DATA_ID_E, 'UniformOutput', false);
    DATA_ID_M = cellfun(@(x) x(1:Truncate), DATA_ID_M, 'UniformOutput', false);
end

IdUnique_E = unique(DATA_ID_E,'stable');
if length(IdUnique_E) < length(DATA_ID_E)
    warning('WARNING! Not all input IDs are unique in DATA_ID_E')
end

IdUnique_M = unique(DATA_ID_M,'stable');
if length(IdUnique_M) < length(DATA_ID_M)
    warning('WARNING! Not all input IDs are unique in DATA_ID_M')
end

[~,ia,ib] = intersect(DATA_ID_E,DATA_ID_M);

DATA_E = EditSamplesDATA(DATA_E,DATA_E.RowId(ia),'keep','stable');
DATA_M = EditSamplesDATA(DATA_M,DATA_M.RowId(ib),'keep','stable');