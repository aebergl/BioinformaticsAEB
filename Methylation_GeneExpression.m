function [DATA DATA_ME]= Methylation_GeneExpression(DATA_M,IdM,CpG_List,DATA_E,IdE,Gene_List)

Truncate = false;

%% Align Datasets
%Get Ids

IdM = ConvertStr(IdM,'string');
IdE = ConvertStr(IdE,'string');

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

Variables=["nTotal CPGs" "min r Spearman" "min r Pearson" "n <0.5 Spearman" "n <0.5 Pearson" "Average r Spearman" "Average r Pearson" "Mean Methylation" "Methylation Range" "Mean Expression" "Expression Range"];
DATA = CreateDataStructure(DATA_E.nCol,length(Variables),[],[]);
DATA.RowId = DATA_E.ColId;
DATA.RowAnnotation = DATA_E.ColAnnotation;
DATA.RowAnnotationFields = DATA_E.ColAnnotationFields;
DATA.ColId=Variables;

GeneNameIndx = 5;
GeneIdColumnNames = {'UCSC_RefGene_Name','genesUniq'};
GeneIdColumnNames = {'genesUniq'};
ChrPosColumnName = 'CpG_beg';

x1 = ones(DATA_E.nCol,1) * NaN;
x2 = ones(DATA_E.nCol,1) * NaN;
x3 = ones(DATA_E.nCol,1) * NaN;
x4 = ones(DATA_E.nCol,1) * NaN;
x5 = ones(DATA_E.nCol,1) * NaN;
x6 = ones(DATA_E.nCol,1) * NaN;
x7 = ones(DATA_E.nCol,1) * NaN;
x8 = ones(DATA_E.nCol,1) * NaN;
x9 = ones(DATA_E.nCol,1) * NaN;
x10 = ones(DATA_E.nCol,1) * NaN;
x11 = ones(DATA_E.nCol,1) * NaN;
r_cutoff = -0.5;
GeneColumns = DATA_M.ColAnnotation(:,35);
GeneDATA = cellfun(@(x) strsplit(x,';'), GeneColumns, 'UniformOutput', 0);

X_E = ones(DATA_E.nRow,DATA_E.nCol) * NaN;
X_M = ones(DATA_E.nRow,DATA_E.nCol) * NaN;
ColId_E = strings(DATA_E.nCol,1);
ColId_M = strings(DATA_E.nCol,1);
parfor i = 1:DATA_E.nCol
    Gene_Id = DATA_E.ColAnnotation(i,GeneNameIndx);
    if ~matches(Gene_Id,["NA" "---" ""])
        x_gene = DATA_E.X(:,i);
        DATA_tmp = FindCpGsForGene(DATA_M,Gene_Id,GeneIdColumnNames,ChrPosColumnName,'PreSplit',GeneDATA);
        if ~isempty(DATA_tmp)
            r_Spearman = corr(x_gene,DATA_tmp.X,'rows','pairwise','Type','Spearman');
            r_Pearson = corr(x_gene,DATA_tmp.X,'rows','pairwise','Type','Pearson');
            indx_Spearman = r_Spearman < r_cutoff;
            indx_Pearson = r_Pearson < r_cutoff;
            x1(i) = length(r_Spearman);
            x2(i) = min(r_Spearman);
            x3(i) = min(r_Pearson);
            x4(i) = sum(indx_Spearman);
            x5(i) = sum(indx_Pearson);
            x_methylation = mean(DATA_tmp.X(:,indx_Spearman),2,'omitnan');
            x6(i) = corr(x_gene,x_methylation,'rows','pairwise','Type','Spearman');
            x7(i) = corr(x_gene,mean(DATA_tmp.X(:,indx_Pearson),2,'omitnan'),'rows','pairwise','Type','Pearson');
            x8(i) = mean(x_methylation,'omitnan');
            x9(i) = range(x_methylation,'omitnan');
            x10(i) = mean(x_gene,'omitnan');
            x11(i) = range(x_gene,'omitnan');
            X_E(:,i)  = x_gene;
            X_M(:,i)  = x_methylation;
            ColId_E(i) = append(Gene_Id, " Expression");
            ColId_M(i) = append(Gene_Id, " Methylation");
        end
    end
end
DATA.X = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11];

DATA_ME = CreateDataStructure(DATA_E.nRow,DATA_E.nCol*2,[],[]);
DATA_ME.RowId = DATA_E.RowId;
DATA_ME.RowAnnotation = DATA_E.RowAnnotation;
DATA_ME.RowAnnotationFields = DATA_E.RowAnnotationFields;
DATA_ME.ColId = [ColId_M;ColId_E];
DATA_ME.X = [X_M X_E];








