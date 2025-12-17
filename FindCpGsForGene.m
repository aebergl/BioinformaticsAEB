function DATA = FindCpGsForGene(DATA,GeneId,GeneIdColumnNames,ChrPosColumnName,varargin)

i=0;
str_cell = [];
while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'PreSplit')
        i = i + 1;
        str_cell = varargin{i};
    end
end

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
    for j=1:length(GeneId)
        GeneIdToUse = GeneId(j);
        if isempty(str_cell)
            indx_col = strcmp(GeneIdColumnNames(i),DATA.ColAnnotationFields);
            str_cell = cellfun(@(x) strsplit(x,';'), DATA.ColAnnotation(:,indx_col), 'UniformOutput', 0);
        end    
        indx_cell = cellfun(@(x) matches(x,GeneIdToUse), str_cell, 'UniformOutput', 0);
        indx = cellfun(@(x) any(x), indx_cell, 'UniformOutput', 1);
        indx_gene = indx_gene | indx;
    end
end
if ~any(indx_gene)
    warning('No CpG probes found %s',join(GeneId))
    DATA = [];
    return
end

DATA =  EditVariablesDATA(DATA,DATA.ColId(indx_gene),'Keep');
ChrPos = DATA.ColAnnotation(:,ChrPosColumn);
ChrPos = cellfun(@(x) str2double(x), ChrPos, 'UniformOutput', 0);
ChrPos = cell2mat(ChrPos);
ChrPos = double(ChrPos);
[~,sortindx] = sort(ChrPos);

DATA =  EditVariablesDATA(DATA,DATA.ColId(sortindx),'Keep','Stable');