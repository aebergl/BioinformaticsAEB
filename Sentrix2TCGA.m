function TCGA_ID = Sentrix2TCGA(SentrixId)

FilePath = '/Users/bergluae/AEBERGL/DATA/METHYLATION/TCGA_NOOB/Sentrix2TCGAId';

c = dir([FilePath,filesep,'*.sdrf.txt']);
FileList = {c.name};
nFiles = numel(FileList);

nVal = length(SentrixId);
TCGA_ID = cell(nVal,1);
TCGA_ID(:) = {'---'};
nFound=0;

for i=1:nFiles
    FileName = FileList{i};
    FullFileName = [FilePath,filesep,FileName];
    opts = detectImportOptions(FullFileName,'FileType','delimitedtext','NumHeaderLines',0,'ReadVariableNames',true,'PreserveVariableNames',true);
    [SelectedVariables]  = intersect(opts.VariableNames,{'Comment [TCGA Barcode]','Array Data File'},'Stable');
    opts.SelectedVariableNames = SelectedVariables;
    C = readcell(FullFileName,opts);
    Sid=C(:,2);
    Tid=C(:,1);
    indx=contains(Sid,'_Red');
    Sid(indx)=[];
    Tid(indx)=[];
    Sid=strrep(Sid,'_Grn.idat','');
    [indx_SentrixId,indx_Sid]  = ismember(SentrixId,Sid);
    indx_Sid = indx_Sid(indx_Sid>0);
    TCGA_ID(indx_SentrixId) = Tid(indx_Sid);
    if sum(strcmp('---',TCGA_ID)) == 0
       break
    end
end