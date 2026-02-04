function DATA_Out = Methylation_Generate_Gene_CpG_Structure(DATA,GeneList,GeneColumnName)

DATA_Out = [];


GeneIdColumnNames = {'UCSC_RefGene_Name','genesUniq'};
ChrPosColumnName  ='CpG_beg';
ChrNameColumnName = 'CpG_chrm';


GeneList = ConvertStr(GeneList,'string');
GeneIdColumnNames = ConvertStr(GeneIdColumnNames,'string');
ChrPosColumnName = ConvertStr(ChrPosColumnName,'string');

GeneList = unique(GeneList);
GeneList(matches(GeneList,["NA" "---" ""])) = [];

ChrPosColumn = strcmpi(ChrPosColumnName,DATA.ColAnnotationFields);
if ~any(ChrPosColumn)
    error('Given value for ChrPos not found, %s not found in ColAnnotationFields',ChrPosColumnName)
end

ChrNameColumn = strcmpi(ChrNameColumnName,DATA.ColAnnotationFields);
if ~any(ChrNameColumn)
    error('Given value for Chromosome not found, %s not found in ColAnnotationFields',ChrNameColumnName)
end


GeneSymbolColumn = matches(DATA.ColAnnotationFields,GeneIdColumnNames);
if ~any(GeneSymbolColumn)
    error('Given value for Gene Symbol not found, %s not found in ColAnnotationFields',GeneIdColumnNames(:))
end

StrCell = cell(length(GeneIdColumnNames),1);
for i=1:length(GeneIdColumnNames)
    indx_col = strcmp(GeneIdColumnNames(i),DATA.ColAnnotationFields);
    StrCell{i} = cellfun(@(x) strsplit(x,';'), DATA.ColAnnotation(:,indx_col), 'UniformOutput', 0);
end

GeneListOut = strings(length(GeneList),1);
CpGProbes = cell(length(GeneList),1);

parfor j=1:length(GeneList)
    GeneIdToUse = GeneList(j);
    indx_gene = false(DATA.nCol,1);
    for i=1:length(GeneIdColumnNames)
        indx_cell = cellfun(@(x) matches(x,GeneIdToUse), StrCell{i}, 'UniformOutput', 0);
        indx = cellfun(@(x) any(x), indx_cell, 'UniformOutput', 1);
        indx_gene = indx_gene | indx;


    end
    if any(indx_gene)
        % Check that all the CpG are on the same Chromosome
        Chr = DATA.ColAnnotation(indx_gene,ChrNameColumn);
        if length(unique(Chr)) == 1
            ChrPos = DATA.ColAnnotation(indx_gene,ChrPosColumn);
            ChrPos = cellfun(@(x) str2double(x), ChrPos, 'UniformOutput', 0);
            ChrPos = cell2mat(ChrPos);
            ChrPos = double(ChrPos);
            [~,sortindx] = sort(ChrPos);
            CpGs = DATA.ColId(indx_gene);
            CpGs = CpGs(sortindx);
            GeneListOut(j) = GeneIdToUse;
            CpGProbes{j} = string(CpGs);
        else
            GeneIdToUse
            Chr
        end

    else
        %fprintf('%s Not found\n',GeneIdToUse)
    end
end
indx = ~matches(GeneListOut,"");
DATA_Out.GeneColumnName = GeneColumnName;
DATA_Out.GeneList = GeneListOut(indx);
DATA_Out.CpGProbes = CpGProbes(indx);

% DATA_Out.GeneList = GeneListOut(1:counter);
% DATA_Out.CpGProbes = CpGProbes(1:counter);