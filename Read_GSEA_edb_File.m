function DATA_out = Read_GSEA_edb_File(InputFile,DATA,ResultName)

if ~isfile(InputFile)
    error('Could not open %s',InputFile)
end

S = readstruct(InputFile,'FileType','xml');

VarName = {'SIZE','ES','NES','NOM p-val','FDR q-val','FWER p-val','RANK AT MAX','RANK SCORE'};
nRow = length(S.DTG);
nCol = length(VarName);

DATA_out = CreateDataStructure(nRow,nCol,[],[]);

C = struct2cell(S.DTG);
DATA_out.X(:,2:end)=squeeze(cell2mat(C([4:8 13:14],1,:)))';

nGenes = squeeze((C(11,1,:)));
nGenes = cellfun(@(x) split(x),nGenes,'UniformOutput',false);
nGenes = cellfun(@(x) length(x),nGenes,'UniformOutput',false);
nGenes = cell2mat(nGenes);
DATA_out.X(:,1) = nGenes;

Names = squeeze((C(3,1,:)));
Names = [Names{:}]';
Names = strrep(Names,'gene_sets.gmt#','');

DATA_out.ColId = string(VarName');
DATA_out.RowId = strcat(ResultName," ",Names);
