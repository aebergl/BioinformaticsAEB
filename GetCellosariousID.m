function [RESULTS, OutputFields] = GetCellosariousID(IdList,CelloDATA,varargin)
% USAGE:
% [RESULTS, OutputFields] = GetCellosariousID(IdList,CelloDATA,varargin)
% Finds the correct identifier for IdList using Cellosaurus data
%
% INPUTS:
% IdList:       List of Cell Libe Ids to be converted
% CelloDATA:    Data structure with Cellosaurus data
%
% OUTPUTS:
% RESULTS: Figure handle to dot plot figure
%
%   options ---------------------------------------
%
%   'OutputFields'      Cell with names ot output fields to export, [] uses all
%   'ManualCorrection'  N*2 with orig Id and New Id to be used instead
%   'HumanOnly'         true/false if only human cell lines should be exported [true]
%   'CleanChar'         Cell with charatcers to be removed when cleaning [{' ','-', '_'}]
%   'PrintResults'      File name for exported results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Anders Berglund, 2025 aebergl at gmail.com                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OutputFields =      [];
ManualCorrection =  [];
HumanOnly =         true;
CleanChar =         {' ','-', '_'};
ReplaceChar =       "";
PrintResults =      false;
% Check input
if nargin > 2
    ArgsList = {'OutputFields','ManualCorrection','HumanOnly','CleanChar','PrintResults'};
    for j=1:2:numel(varargin)

        ArgType = varargin{j};
        ArgVal = varargin{j+1};
        if ~strncmpi(ArgType,ArgsList,numel(ArgType))
            error('Invalid input option: %s', ArgType);
        else
            switch lower(ArgType)
                case 'outputfields'
                    OutputFields = ArgVal;
                case 'manualcorrection'
                    ManualCorrection = ArgVal;
                case 'humanonly'
                    HumanOnly =ArgVal;
                case 'cleanchar'
                    CleanChar = ArgVal;
                case 'printresults'
                    PrintResults = true;
                    FileOut = ArgVal;
            end

        end
    end
end

n = length(IdList);

if isempty(OutputFields)
    OutputFields = CelloDATA.ColumnFields;
end

[~, OutputFieldsIndx, ResultsIndx] = intersect(OutputFields,CelloDATA.ColumnFields,'stable');

if length(OutputFields) ~= length(OutputFieldsIndx)
    warning('Not all OutputFields found in CELLOSAURUS.ColumnFields')
end

OutputFields(end+1) = {'Match Type'};

OutputFields = [{'Orig Id'},OutputFields];

RESULTS = strings(n,length(OutputFields));


PrimaryId = CelloDATA.CellData(:,strcmp('ID',CelloDATA.ColumnFields));
PrimaryIdClean = replace(PrimaryId,CleanChar,ReplaceChar);

% Start checking

if HumanOnly
    indx_human = contains(CelloDATA.CellData(:,strcmp('OX',CelloDATA.ColumnFields)),'NCBI_TaxID=9606');
else
    indx_human = true(CelloDATA.nCellLines,1);
end
EntryNumHuman = find(indx_human);

for i = 1:n
    Id = IdList(i);
    Id_Clean = replace(Id,CleanChar,ReplaceChar);
    RESULTS(i,1) = Id;

    % Check if Id is in Manual Correction array
    if ~isempty(ManualCorrection)
        indx_MC = strcmp(Id,ManualCorrection(:,1));
        if any(indx_MC)
            Id = ManualCorrection(indx_MC,2);
        end
    end

    % Start to check against Primary Id, case sensative
    indx = strcmp(Id,PrimaryId) & indx_human;

    if any(indx) % Should always be une unique hit
        RESULTS(i,2:end-1) = CelloDATA.CellData(indx,ResultsIndx);
        RESULTS(i,end) = "Primary Id match";
    else

        % Check agains Synonyms, case sensative
        indx = strcmp(Id,CelloDATA.Synonyms);
        EntryNum = CelloDATA.SynonymsIndx(indx);
        if (sum(indx) == 1 || isscalar(unique((CelloDATA.SynonymsIndx(indx))))) && ismember(EntryNum,EntryNumHuman)% Only one unique Hit
            RESULTS(i,2:end-1) = CelloDATA.CellData(EntryNum,ResultsIndx);
            RESULTS(i,end) = "Unique match to synonym";

        elseif sum(indx) > 1 % Multiple hits
            % Check if only one is human
            OX = CelloDATA.CellData(EntryNum,strcmp('OX',CelloDATA.ColumnFields));
            indx_hum = strcmp('NCBI_TaxID=9606; ! Homo sapiens (Human)',OX);
            if sum(indx_hum) == 1
                RESULTS(i,2:end-1) = CelloDATA.CellData(EntryNum(indx_hum),ResultsIndx);
                RESULTS(i,end) = "Multiple match to synonym, only one human";

            else
                indx_id = strcmp(Id_Clean,PrimaryIdClean(EntryNum));
                if sum(indx_id) == 1
                    RESULTS(i,2:end-1) = CelloDATA.CellData(EntryNum(indx_id),ResultsIndx);
                    RESULTS(i,end) = "Multiple match to synonym, one match cleaned primary Id";
                else
                    RESULTS(i,end) = "MULTIPLE MATCHES!";

                end

            end
        else
            indx = strcmpi(Id,PrimaryId) & indx_human;
            if sum(indx) == 1 % Only one unique Hit
                RESULTS(i,2:end-1) = CelloDATA.CellData(indx,ResultsIndx);
                RESULTS(i,end) = "Unique match to Primary Id, ignoring case";
            else
                indx = strcmpi(Id,CelloDATA.Synonyms); % Case insensative
                if sum(indx) == 1 || isscalar(unique((CelloDATA.SynonymsIndx(indx))))% Only one unique Hit
                    EntryNum = CelloDATA.SynonymsIndx(indx);
                    RESULTS(i,2:end-1) = CelloDATA.CellData(EntryNum(1),ResultsIndx);
                    RESULTS(i,end) = "Unique match to synonym, ignoring case";
                else
                    indx = strcmpi(Id_Clean,PrimaryIdClean) & indx_human ;
                    if sum(indx) == 1 % Only one unique Hit
                        RESULTS(i,2:end-1) = CelloDATA.CellData(indx,ResultsIndx);
                        RESULTS(i,end) = "Primary Id matched with cleaned Id";
                    else
                        RESULTS(i,end) = "NO MATCH FOUND";
                    end
                end
            end
        end
    end
end

if PrintResults
    [fid,message] = fopen(FileOut,'w');
    if  fid == -1
        disp(FileOut)
        disp(message)
        return
    end
    fprintf(fid,'%s\t',OutputFields{:});
    fprintf(fid,'\n');
    for i = 1:n
        fprintf(fid,'%s\t',RESULTS(i,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
end
