function RESULTS = Methylation_Calc_Methylation_levels(DATA,CHR_Column)

ChrToExclude = ["chrX" "chrY" "NA"];

VarNames = {'Average Methylation','Ratio >0.5/<0.5','Ratio >0.7/<0.3','Ratio >0.75/<0.25','Ratio >0.8/<0.2'};
nVar = length(VarNames);


indx_column = strcmpi(CHR_Column,DATA.ColAnnotationFields);
if ~any(indx_column)
    error('Error. \n%s not found in DATA.ColAnnotationFields',CHR_Column);
elseif sum(indx_column) > 1
    error('Warning. \nMultiple matches for %s found in DATA.ColAnnotationFields',CHR_Column);
end


% only use cg_probes
indx_include = strncmpi("cg*",DATA.ColId,2);

% Get indx for CpG-probes to remove
for i=1:length(ChrToExclude)
    indx_tmp = strcmpi(ChrToExclude(i),DATA.ColAnnotation(:,indx_column));
    indx_include(indx_tmp) = false;
end

% Create results DATA structure
RESULTS = CreateDataStructure(DATA.nRow,nVar,[],[]);

RESULTS.Title = 'Methylation Level';
RESULTS.Info.Source = inputname(1);

RESULTS.ColId=VarNames;
RESULTS.RowId = DATA.RowId;
RESULTS.RowAnnotationFields = DATA.RowAnnotationFields;
RESULTS.RowAnnotation = DATA.RowAnnotation;
if isfield(DATA,'SURVIVAL')
    RESULTS.SURVIVAL = DATA.SURVIVAL;
end



% Average Methylation
RESULTS.X(:,1) = mean(DATA.X(:,indx_include),2,'omitnan');

% Ratio >0.5/<0.5
RESULTS.X(:,2) = (sum(DATA.X(:,indx_include)>0.5,2) - sum(DATA.X(:,indx_include)<0.5,2)) ./ (sum(DATA.X(:,indx_include)>0.5,2) + sum(DATA.X(:,indx_include)<0.5,2));

% Ratio >0.7/<0.3
%RESULTS.X(:,3) = sum(DATA.X(:,indx_include)>0.7,2) ./ sum(DATA.X(:,indx_include)<0.3,2);
RESULTS.X(:,3) = (sum(DATA.X(:,indx_include)>0.7,2) - sum(DATA.X(:,indx_include)<0.3,2)) ./ (sum(DATA.X(:,indx_include)>0.7,2) + sum(DATA.X(:,indx_include)<0.3,2));

% Ratio >0.75/<0.25
%RESULTS.X(:,4) = sum(DATA.X(:,indx_include)>0.75,2) ./ sum(DATA.X(:,indx_include)<0.25,2);
RESULTS.X(:,4) = (sum(DATA.X(:,indx_include)>0.75,2) - sum(DATA.X(:,indx_include)<0.25,2)) ./ (sum(DATA.X(:,indx_include)>0.75,2) + sum(DATA.X(:,indx_include)<0.25,2));

% Ratio >0.8/<0.2
%RESULTS.X(:,5) = sum(DATA.X(:,indx_include)>0.8,2) ./ sum(DATA.X(:,indx_include)<0.2,2);
RESULTS.X(:,5) = (sum(DATA.X(:,indx_include)>0.8,2) - sum(DATA.X(:,indx_include)<0.2,2)) ./ (sum(DATA.X(:,indx_include)>0.8,2) + sum(DATA.X(:,indx_include)<0.2,2));
