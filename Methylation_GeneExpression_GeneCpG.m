function [DATA, DATA_ME, DATA_ALL]= Methylation_GeneExpression_GeneCpG(DATA_M,IdM,GeneCpG,DATA_E,IdE)

Truncate = false;

%% Align Datasets
%Get Ids

IdM = ConvertStr(IdM,'string');
IdE = ConvertStr(IdE,'string');

GeneSymbolColumn = strcmp(GeneCpG.GeneIdColumn,DATA_E.ColAnnotationFields);
if ~any(GeneSymbolColumn)
    error('Given value for Gene Symbol not found, %s not found in ColAnnotationFields',GeneCpG.GeneIdColumn)
end


if isempty(IdM)
    M_Id = DATA_M.RowId;
else
    indx_IdM = strcmpi(IdM,DATA_M.RowAnnotationFields);
    if ~any(indx_IdM)
        error('Error. \n%s not found in DATA_M.RowAnnotationFields',IdM);
    elseif sum(indx_IdM) > 1
        error('Warning. \nMultiple matches for %s found in DATA_M.RowAnnotationFields',IdM);
    else
        M_Id = DATA_M.RowAnnotation(:,indx_IdM);
    end

end

if isempty(IdE)
    E_Id = DATA_E.RowId;
else
    indx_IdE = strcmpi(IdE,DATA_E.RowAnnotationFields);
    if ~any(indx_IdE)
        error('Error. \n%s not found in DATA_E.RowAnnotationFields',IdE);
    elseif sum(indx_IdE) > 1
        error('Warning. \nMultiple matches for %s found in DATA_E.RowAnnotationFields',IdE);
    else
        E_Id = DATA_E.RowAnnotation(:,indx_IdE);
    end

end

if Truncate
    File_Id = cellfun(@(x) x(1:Truncate), File_Id, 'UniformOutput', false);
    DATA_Id = cellfun(@(x) x(1:Truncate), DATA_Id, 'UniformOutput', false);
end
[~,M_indx,E_indx] = intersect(M_Id,E_Id,'stable');

DATA_M = EditSamplesDATA(DATA_M,DATA_M.RowId(M_indx),'Keep','Stable');
DATA_E = EditSamplesDATA(DATA_E,DATA_E.RowId(E_indx),'Keep','Stable');

nGenes = length(GeneCpG.GeneList);

Variables=["nTotal CPGs" "min r Spearman" "min r Pearson" "n <0.5 Spearman" "n <0.5 Pearson" "Average r Spearman" "p Spearman" "Average r Pearson" "p Pearson" "Mean Methylation" "Methylation Range" "Mean Expression" "Expression Range"];
DATA = CreateDataStructure(nGenes,length(Variables),[],[]);
DATA.RowAnnotationFields = DATA_E.ColAnnotationFields;
DATA.ColId=Variables;


RowId = strings(nGenes,1);
RowAnnotation = strings(nGenes, size(DATA_E.ColAnnotation,2));
x1 = ones(nGenes,1) * NaN;
x2 = ones(nGenes,1) * NaN;
x3 = ones(nGenes,1) * NaN;
x4 = ones(nGenes,1) * NaN;
x5 = ones(nGenes,1) * NaN;
x6 = ones(nGenes,1) * NaN;
x7 = ones(nGenes,1) * NaN;
x8 = ones(nGenes,1) * NaN;
x9 = ones(nGenes,1) * NaN;
x10 = ones(nGenes,1) * NaN;
x11 = ones(nGenes,1) * NaN;
x12 = ones(nGenes,1) * NaN;
x13 = ones(nGenes,1) * NaN;

r_cutoff = -0.5;

X_E = ones(DATA_E.nRow,nGenes) * NaN;
X_M = ones(DATA_E.nRow,nGenes) * NaN;
ColId_E = strings(nGenes,1);
ColId_M = strings(nGenes,1);
AllCorrSpearman = cell(nGenes,1);
AllMethylationRange = cell(nGenes,1);

parfor i = 1:nGenes
    Gene_Id = GeneCpG.GeneList(i);
    gene_indx = strcmp(Gene_Id,DATA_E.ColAnnotation(:,GeneSymbolColumn));
    if any(gene_indx)
        x_gene = DATA_E.X(:,gene_indx);
        DATA_tmp = EditVariablesDATA(DATA_M,GeneCpG.CpGProbes{i},'keep','stable');
        if ~isempty(DATA_tmp)
            RowId(i) = Gene_Id;
            RowAnnotation(i,:) = DATA_E.ColAnnotation(gene_indx,:);
            x_methylation = DATA_tmp.X;
            r_Spearman = corr(x_gene,x_methylation,'rows','pairwise','Type','Spearman');
            AllCorrSpearman{i}  = r_Spearman;
            AllMethylationRange{i}  = range(x_methylation,1);
            r_Pearson = corr(x_gene,x_methylation,'rows','pairwise','Type','Pearson');
            indx_Spearman = r_Spearman < r_cutoff;
            indx_Pearson = r_Pearson < r_cutoff;
            x1(i) = length(r_Spearman);
            x2(i) = min(r_Spearman);
            x3(i) = min(r_Pearson);
            x4(i) = sum(indx_Spearman);
            x5(i) = sum(indx_Pearson);
            x_mean_methylation = mean(x_methylation(:,indx_Spearman),2,'omitnan');
            [x6(i), x7(i)] = corr(x_gene,x_mean_methylation,'rows','pairwise','Type','Spearman');
            [x8(i), x9(i)] = corr(x_gene,mean(DATA_tmp.X(:,indx_Pearson),2,'omitnan'),'rows','pairwise','Type','Pearson');
            x10(i) = mean(x_mean_methylation,'omitnan');
            x11(i) = range(x_mean_methylation,'omitnan');
            x12(i) = mean(x_gene,'omitnan');
            x13(i) = range(x_gene,'omitnan');
            X_E(:,i)  = x_gene;
            X_M(:,i)  = x_mean_methylation;
            ColId_E(i) = append(Gene_Id, " Expression");
            ColId_M(i) = append(Gene_Id, " Methylation");
        end
    end
end
X = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13];
indx = ~matches(RowId,"");
DATA.RowId = RowId(indx);
DATA.RowAnnotation = RowAnnotation(indx,:);
DATA.X = X(indx,:);
DATA.nRow = length(DATA.RowId);

DATA_ALL.AllCorrSpearman = AllCorrSpearman(indx);
DATA_ALL.AllMethylationRange = AllMethylationRange(indx);

DATA_ME = CreateDataStructure(DATA_E.nRow,DATA.nRow*2,[],[]);
DATA_ME.RowId = DATA_E.RowId;
DATA_ME.RowAnnotation = DATA_E.RowAnnotation;
DATA_ME.RowAnnotationFields = DATA_E.RowAnnotationFields;
DATA_ME.ColId = [ColId_M(indx);ColId_E(indx)];
DATA_ME.X = [X_M(:,indx) X_E(:,indx)];








