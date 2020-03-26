function ReadTCGA_Survival_File

SurvivalTypes = {'OS','PFI','DSS','DFI'};

C = readcell('TCGA-CDR-SupplementalTableS1.xlsx','Sheet','TCGA-CDR');

%Remove first column with just numbers
C(:,1) = [];

ColumnNames = C(1,:);
C(1,:) = [];

[~,EventIndx] = ismember(SurvivalTypes,ColumnNames);
SurvEvent = C(:,EventIndx);

TimeIndx = EventIndx + 1;

SurvTime = C(:,TimeIndx);
indx_missing = ~cellfun(@(x) isnumeric(x),SurvTime);
[SurvTime{indx_missing}]  = deal(NaN);
SurvTime = cell2mat(SurvTime);


