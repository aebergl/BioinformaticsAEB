function GSEA = Read_GSEA_File(FileName)

% Open GSEA Text file File
[fid,message] = fopen(FileName,'r');
if  fid == -1
    disp(FileName)
    disp(message)
    return
end

fgetl(fid);
fgetl(fid);
fgetl(fid);

fline = fgetl(fid);
tmp=split(fline);
GSEA.Collection = tmp(2);
tmp = textscan(fid,'%s %d',1,'delimiter',{':'});
GSEA.INFO.numPathways = tmp{2};

tmp=textscan(fid,'%s %d',1,'delimiter',{':'});
GSEA.INFO.numGeneSetsInCollection = tmp{2};

tmp=textscan(fid,'%s %d',1,'delimiter',{':'});
GSEA.INFO.numGenesInComparison = tmp{2};

tmp=textscan(fid,'%s %d',1,'delimiter',{':'});
GSEA.INFO.numGenesInUniverse = tmp{2};

fgetl(fid);

fgetl(fid);
fgetl(fid);
tmp=textscan(fid,'%s %d %s %d %f %f %f',GSEA.INFO.numPathways,'delimiter',{'\t'});

GSEA.PATHWAYS.Name = tmp{1};
GSEA.PATHWAYS.numGenesInSet = tmp{2};
GSEA.PATHWAYS.Description = tmp{3};
GSEA.PATHWAYS.numGenesInOveralap = tmp{4};
GSEA.PATHWAYS.Ratio = tmp{5};
GSEA.PATHWAYS.p = tmp{6};
GSEA.PATHWAYS.q = tmp{7};

fgetl(fid);
fgetl(fid);

fgetl(fid);

fgetl(fid);
fgetl(fid);
fgetl(fid);

PrintFormat = strcat('%d%s%s',repmat('%s',1,GSEA.INFO.numPathways));
tmp=textscan(fid,PrintFormat,GSEA.INFO.numGenesInComparison,'delimiter',{'\t'});


GSEA.MATRIX.EntrezId = tmp{1};
GSEA.MATRIX.GeneSymbol = tmp{2};
GSEA.MATRIX.GeneDescription = tmp{3};
GSEA.MATRIX.PathwayNames = GSEA.PATHWAYS.Name;
Mtmp =  tmp(4:end);
GSEA.MATRIX.PathwayMatrix = [Mtmp{:}];


fclose(fid);
