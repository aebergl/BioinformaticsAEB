function DATA = FindCpGsForGene(DATA,GeneId,GeneIdColumnNames,ChrPosColumnName,varargin)

CpGProbesOnly = false;

GeneId = ConvertStr(GeneId,'string');
GeneIdColumnNames = ConvertStr(GeneIdColumnNames,'string');
ChrPosColumnName = ConvertStr(ChrPosColumnName,'string');


ChrPosColumn = strcmpi(ChrPosColumnName,DATA.ColAnnotationFields);
if ~any(ChrPosColumn)
    error('Given value for ChrPos not found, %s not found in ColAnnotationFields',ChrPosColumnName)
end


indx_gene = false(DATA.nCol,1);
for i=1:length(GeneIdColumnNames)
    if length(GeneId) > 1
        GeneIdToUse = GeneId(i);
    else
        GeneIdToUse = GeneId;
    end
    indx_col = strcmp(GeneIdColumnNames(i),DATA.ColAnnotationFields);
    indx_1 = strcmp(GeneIdToUse,DATA.ColAnnotation(:,indx_col));
    indx_2 = contains(DATA.ColAnnotation(:,indx_col),strcat(";",GeneIdToUse));
    indx_3 = contains(DATA.ColAnnotation(:,indx_col),strcat(GeneIdToUse,";"));
    indx = indx_1 | indx_2 | indx_3;
    indx_gene = indx_gene | indx;
end
if ~any(indx_gene)
    error('No CpG probes found %s',join(GeneId))
end

DATA =  EditVariablesDATA(DATA,DATA.ColId(indx_gene),'Keep');
ChrPos = DATA.ColAnnotation(:,ChrPosColumn);
ChrPos = cellfun(@(x) str2num(x), ChrPos, 'UniformOutput', 0);
ChrPos = cell2mat(ChrPos);
ChrPos = double(ChrPos);
[~,sortindx] = sort(ChrPos);

DATA =  EditVariablesDATA(DATA,DATA.ColId(sortindx),'Keep','Stable');