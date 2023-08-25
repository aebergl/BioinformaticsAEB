function LIST = MethylationGetImprintedGenes(SourceTable,AnnotationFile)


%%
% Read Source table
opts_1 = detectImportOptions(SourceTable,'FileType','delimitedtext','TreatAsMissing',{'NA','TBD'});
TableData = readcell(SourceTable,opts_1,"TextType","char");

% Get Imprinted Ids
pat = lettersPattern + "_" + digitsPattern;
IDs = string(TableData(:,1));
IDs = extract(IDs,pat);

nIDs = length(IDs);

% Get Chromosome number and range
GenomicCoords = string(TableData(:,2));
ChrNum = extractBefore(GenomicCoords,':');
StartPos = str2double(extractBetween(GenomicCoords,':','-'));
StopPos = str2double(extractAfter(GenomicCoords,'-'));

% get Parental Origin of Methylation
ParentalOrigin = string(TableData(:,3));

% get Nearest Transcript
NearestTranscript = string(TableData(:,4));


% get Nearest Transcript
DistanceToNT = str2double(string(TableData(:,5)));



%%
% Read Annotation file

opts_2 = detectImportOptions(AnnotationFile,'FileType','delimitedtext','TreatAsMissing',{'NA','TBD'});

CellData = readcell(AnnotationFile,opts_2,"TextType","char");

% Chip Chromosome number
ChipChrNum=string(CellData(:,1));

% Chip Chromosome coordinate
ChipPos=[CellData(:,2)];
indx_missing = ~cellfun(@(x) isnumeric(x),ChipPos);
[ChipPos{indx_missing}]  = deal(NaN);
ChipPos = cell2mat(ChipPos);

% Chip probe Id
ChipProbeID=string(CellData(:,5));

% Chip genesUniq
ChipGenesUniq=string(CellData(:,6));

% Chip TranscriptTypes
ChipTranscriptTypes=string(CellData(:,6));

%%
fprintf('%s\t',opts_1.VariableNames{:})
fprintf('probe Id\tGenes Uniq\tTranscript Type\n')
LIST = "";
for i=1:nIDs
    indx=find((ChipChrNum==ChrNum(i)) & ChipPos > StartPos(i) & ChipPos < StopPos(i));
    if isempty(indx)
            fprintf('%s\t%s\t%s\t%s\t%i\n',IDs(i),GenomicCoords(i),ParentalOrigin(i),NearestTranscript(i),DistanceToNT(i));

    else
        for j=1:length(indx)
            fprintf('%s\t%s\t%s\t%s\t%i\t',IDs(i),GenomicCoords(i),ParentalOrigin(i),NearestTranscript(i),DistanceToNT(i));
            fprintf('%s\t%s\t%s\n',CellData{indx(j),[5 6 8 ]})
            LIST = [LIST;CellData(indx(j),5)];
        end
    end

end

