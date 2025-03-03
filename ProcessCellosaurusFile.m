function DATA = ProcessCellosaurusFile(FileIn) 
Terminator = '//';

%Read already downloaded 'EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv'
if ~isfile(FileIn) || isempty(FileIn)
    try
        warning('Could not find cellosaurus.txt in this directory')
        warning('Starting to download it now, will probably take a couple of minutes')
        websave('cellosaurus.txt','https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt');
        
    catch
        error('Could not load EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv from https://gdc.cancer.gov/about-data/publications/pancanatlas')
    end
    L = readlines('cellosaurus.txt');
else
    L = readlines(FileIn);

end
fprintf('Processing cellosaurus_test.txt\n');
nRow = length(L);

% Get number of Cell lines
nCellLines=sum(strcmp(L,Terminator));

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

S = repmat(struct( 'ID',"",'AC',"",'AS',"",'SY',"",'DR',"",'RX',"",'WW',"",'CC',"",'ST',"",'DI',"",'OX',"",'HI',"",'OI',"",'SX',"",'AG',"",'CA',"",'DT',""),nCellLines,1);

FN = fieldnames(S);
CellCounter = 0;
RowCounter = 0;
while RowCounter < nRow;
    RowCounter = RowCounter + 1;
    tmp = L(RowCounter);
    if strncmp(tmp,'ID',2) %Starts an Entry
        CellCounter = CellCounter + 1;
        S(CellCounter).ID = extractAfter(tmp,5);
        while ~strncmp(tmp,Terminator,2)
            RowCounter = RowCounter + 1;
            tmp = L(RowCounter);
            switch extractBefore(tmp,3)
                case 'AC'
                    S(CellCounter).AC = extractAfter(tmp,5);
                case 'AS'
                    S(CellCounter).AS = extractAfter(tmp,5);
                case 'SY'
                    S(CellCounter).SY = extractAfter(tmp,5);
                case 'DR'
                    if strlength(S(CellCounter).DR)
                        S(CellCounter).DR = strcat(S(CellCounter).DR,"|",extractAfter(tmp,5));
                    else
                        S(CellCounter).DR = extractAfter(tmp,5);
                    end
                case 'RX'
                    if strlength(S(CellCounter).RX)
                        S(CellCounter).RX = strcat(S(CellCounter).RX,"|",extractAfter(tmp,5));
                    else
                        S(CellCounter).RX = extractAfter(tmp,5);
                    end
                case 'WW'
                    if strlength(S(CellCounter).WW)
                        S(CellCounter).WW = strcat(S(CellCounter).WW,"|",extractAfter(tmp,5));
                    else
                        S(CellCounter).WW = extractAfter(tmp,5);
                    end
                case 'CC'
                    if strlength(S(CellCounter).CC)
                        S(CellCounter).CC = strcat(S(CellCounter).CC,"|",extractAfter(tmp,5));
                    else
                        S(CellCounter).CC = extractAfter(tmp,5);
                    end
                case 'ST'
                    if strlength(S(CellCounter).ST)
                        S(CellCounter).ST = strcat(S(CellCounter).ST,"|",extractAfter(tmp,5));
                    else
                        S(CellCounter).ST = extractAfter(tmp,5);
                    end
                case 'DI'
                    if strlength(S(CellCounter).DI)
                        S(CellCounter).DI = strcat(S(CellCounter).DI,"|",extractAfter(tmp,5));
                    else
                        S(CellCounter).DI = extractAfter(tmp,5);
                    end
                case 'OX'
                    if strlength(S(CellCounter).OX)
                        S(CellCounter).OX = strcat(S(CellCounter).OX,"|",extractAfter(tmp,5));
                    else
                        S(CellCounter).OX = extractAfter(tmp,5);
                    end
                case 'HI'
                    if strlength(S(CellCounter).HI)
                        S(CellCounter).HI = strcat(S(CellCounter).HI,"|",extractAfter(tmp,5));
                    else
                        S(CellCounter).HI = extractAfter(tmp,5);
                    end
                case 'OI'
                    if strlength(S(CellCounter).OI)
                        S(CellCounter).OI = strcat(S(CellCounter).OI,"|",extractAfter(tmp,5));
                    else
                        S(CellCounter).OI = extractAfter(tmp,5);
                    end
                case 'SX'
                    S(CellCounter).SX = extractAfter(tmp,5);
                case 'AG'
                    S(CellCounter).AG = extractAfter(tmp,5);
                case 'CA'
                    S(CellCounter).CA = extractAfter(tmp,5);
                case 'DT'
                    S(CellCounter).DT = extractAfter(tmp,5);
            end


        end

    end


end
DATA.CellLine=S;
DATA.nCellLines = nCellLines;

% Process Synonyms for easy search
SynonymsRaw = [S(:).SY]';
SynonymsSplit = arrayfun(@(x) strsplit(x,'; '), SynonymsRaw, 'UniformOutput', 0);
% Get indx CELL
Indx_Num = cellfun(@(x) strlength(x)>0, SynonymsSplit, 'UniformOutput', 0);
Indx_Num=cellfun(@(x,y) x.*y, Indx_Num, num2cell([1:length(SynonymsRaw)]'), 'UniformOutput', 0);

% Expand to single vector
Synonyms = [SynonymsSplit{:}]';
Indx_Num = [Indx_Num{:}]';

Synonyms = Synonyms(Indx_Num>0);
Indx_Num = Indx_Num(Indx_Num>0);
DATA.Synonyms = Synonyms;
DATA.SynonymsIndx = Indx_Num;


