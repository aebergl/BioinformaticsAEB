function DATA = ProcessCellosaurusFile(FileIn)
Terminator = "//";
DATA.DateDownloaded = "";
%Read already downloaded 'EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv'
if isempty(FileIn) || ~isfile(FileIn)
    try
        warning('Could not find cellosaurus.txt in this directory')
        warning('Starting to download it now, will probably take a couple of minutes')
        websave('cellosaurus.txt','https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt');

    catch
        error('Could not load cellosaurus.txt from https://ftp.expasy.org/databases/cellosaurus/save ')
    end
    L = readlines('cellosaurus.txt');
    DATA.DateDownloaded = string(datetime('today'));
else
    L = readlines(FileIn);
end
DATA.Version = L(7);
DATA.LastUpdate = L(8);


% ---------  ------------------------------  -----------------------
% Line code  Content                         Occurrence in an entry
% ---------  ------------------------------  -----------------------
% ID         Identifier (cell line name)     Once; starts an entry
% AC         Accession (CVCL_xxxx)           Once
% AS         Secondary accession number(s)   Optional; once
% SY         Synonyms                        Optional; once
% DR         Cross-references                Optional; once or more
% RX         References identifiers          Optional: once or more
% WW         Web pages                       Optional; once or more
% CC         Comments                        Optional; once or more
% ST         STR profile data                Optional; twice or more
% DI         Diseases                        Optional; once or more
% OX         Species of origin               Once or more
% HI         Hierarchy                       Optional; once or more
% OI         Originate from same individual  Optional; once or more
% SX         Sex of cell                     Optional; once
% AG         Age of donor at sampling        Optional; once
% CA         Category                        Once
% DT         Date (entry history)            Once
% //         Terminator                      Once; ends an entry


fprintf('Processing cellosaurus.txt\n');
nRow = length(L);

% Get number of Cell lines
nCellLines=sum(strcmp(L,Terminator));


% Create string array to store all data
ColumnFields = ["ID","AC","AS","SY","DR","RX","WW","CC","ST","DI","OX","HI","OI","SX","AG","CA","DT"];
DATA.ColumnFields = ColumnFields;
DATA.CellData = strings(nCellLines,length(ColumnFields));


CellCounter = 0;
RowCounter = 0;
while RowCounter < nRow
    RowCounter = RowCounter + 1;
    tmp = L(RowCounter);
    if strncmp(tmp,'ID',2) %Starts an Entry
        CellCounter = CellCounter + 1;
        while ~strncmp(tmp,Terminator,2)
            extractBefore(tmp,3);
            pos = strcmp(extractBefore(tmp,3),ColumnFields);
            if strlength(DATA.CellData(CellCounter,pos))
                DATA.CellData(CellCounter,pos) = strcat(DATA.CellData(CellCounter,pos),"|",extractAfter(tmp,5));
            else
                DATA.CellData(CellCounter,pos) = extractAfter(tmp,5);
            end
            RowCounter = RowCounter + 1;
            tmp = L(RowCounter);
        end
    end
end

DATA.nCellLines = nCellLines;

% Process Synonyms for easy search

SynonymsRaw = DATA.CellData(:,4);
SynonymsSplit = arrayfun(@(x) strsplit(x,'; '), SynonymsRaw, 'UniformOutput', 0);

% Get indx CELL
Indx_Num = cellfun(@(x) strlength(x)>0, SynonymsSplit, 'UniformOutput', 0);
Indx_Num = cellfun(@(x,y) x.*y, Indx_Num, num2cell((1:length(SynonymsRaw))'), 'UniformOutput', 0);

% Expand to single vector
Synonyms = [SynonymsSplit{:}]';
Indx_Num = [Indx_Num{:}]';

Synonyms = Synonyms(Indx_Num>0);
Indx_Num = Indx_Num(Indx_Num>0);
DATA.Synonyms = Synonyms;
DATA.SynonymsIndx = Indx_Num;


