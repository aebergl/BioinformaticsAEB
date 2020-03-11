function DATA = AddREFseqGeneInfo(DATA,FileName,ReplaceAppend)

if isempty(FileName)
    ftpobj = ftp('ftp.ncbi.nlm.nih.gov');
    cd(ftpobj,'refseq/H_sapiens/');
    mget(ftpobj,'Homo_sapiens.gene_info.gz');
    gunzip('Homo_sapiens.gene_info.gz');
    [~,file_name,~] = fileparts('Homo_sapiens.gene_info');
    [FidInputFile,message] = fopen('Homo_sapiens.gene_info','r');
    if  FidInputFile == -1
        disp(InputFile)
        disp(message)
        return
    end
    
else
    [~,file_name,~] = fileparts(FileName);
    [FidInputFile,message] = fopen(FileName,'r');
    if  FidInputFile == -1
        disp(FileName)
        disp(message)
        return
    end
    
    
end
tline = fgetl(FidInputFile);
tline =  textscan(tline,'%s','delimiter','\t');
ColumnIds = tline{1};
numColumns = numel(ColumnIds);
Annotation = cell(DATA.nCol,8);
Annotation(:) = {'---'};

S = textscan(FidInputFile,repmat('%s',1,numColumns),'delimiter','\t');

GeneId = S{2};
GeneColumnsToUse = [2 3 9 5 7 8 10 6];
GeneInfo = cat(2,S{GeneColumnsToUse});

[~,indx1,indx2] = intersect(DATA.ColId,GeneId,'Stable');

whos

if length(indx1) < DATA.nCol
    fprintf('WARNING!!! %u columns are missing annotation\n',DATA.nCol - length(indx1));  
end

Annotation(indx1,:) = GeneInfo(indx2,:);
Annotation(cellfun('isempty',Annotation)) = {''};
switch lower(ReplaceAppend)
    case 'replace'
        DATA.ColAnnotationFields = ColumnIds(GeneColumnsToUse);
        DATA.ColAnnotation = Annotation;
    case 'append'
        DATA.ColAnnotationFields = [DATA.ColAnnotationFields; ColumnIds(GeneColumnsToUse)];
        DATA.ColAnnotation = [DATA.ColAnnotation Annotation];   
end















