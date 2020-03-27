function SURVIVAL= ReadTCGA_Survival_File

% OS: overall survival event, 1 for death from any cause, 0 for alive.
% OS.time: overall survival time in days, last_contact_days_to or death_days_to, whichever is larger.

% DSS: disease-specific survival event, 1 for patient whose vital_status was Dead and tumor_status was WITH TUMOR. If a patient died from the disease shown in field of cause_of_death, the status of DSS would be 1 for the patient.  0 for patient whose vital_status was Alive or whose vital_status was Dead and tumor_status was  TUMOR FREE. This is not a 100% accurate definition but is the best we could do with this dataset. Technically a patient could be with tumor but died of a car accident and therefore incorrectly considered as an event.
% DSS.time: disease-specific survival time in days, last_contact_days_to or death_days_to, whichever is larger.
% 
% DFI: disease-free interval event, 1 for patient having new tumor event whether it is a local recurrence, distant metastasis, new primary tumor of the cancer, including cases with a new tumor event whose type is N/A.  Disease free was defined by: first, treatment_outcome_first_course is "Complete Remission/Response"; if the tumor type doesn't have "treatment_outcome_first_course" then disease-free was defined by the value "R0" in the field of "residual_tumor"; otherwise, disease-free was defined by the value "negative" in the field of "margin_status". If the tumor type did not have any of these fields, then its DFI was NA.
% DFI.time: disease-free interval time in days, new_tumor_event_dx_days_to for events, or for censored cases, either last_contact_days_to or death_days_to, whichever is applicable.
%
% PFI: progression-free interval event, 1 for patient having new tumor event whether it was a progression of disease, local recurrence, distant metastasis, new primary tumors all sites , or died with the cancer without new tumor event, including cases with a new tumor event whose type is N/A.
% PFI.time: progression-free interval time in days, for events, either new_tumor_event_dx_days_to or death_days_to,  whichever is applicable; or for censored cases, either last_contact_days_to or death_days_to, whichever is applicable.
% 
C = readcell('TCGA-CDR-SupplementalTableS1.xlsx','Sheet','TCGA-CDR');

%Remove first column with just numbers
C(:,1) = [];

% Collect Column names and remove the first rwo
ColumnNames = C(1,:);
C(1,:) = [];

% Create Survival Structure

SURVIVAL.bcr_patient_barcode = C(:,1);
SURVIVAL.SurvivalTypes = {'OS','PFI','DSS','DFI'};
% Get the event columns
[~,EventIndx] = ismember(SURVIVAL.SurvivalTypes,ColumnNames);
SurvEvent_RAW = C(:,EventIndx);
SURVIVAL.SurvEvent = cell(size(SurvEvent_RAW));
SURVIVAL.SurvEvent(:) = {'NA'};

% Convert OS events
SURVIVAL.SurvEvent(cellfun(@(x) x==1, SurvEvent_RAW(:,1)),1) = {'Dead'};
SURVIVAL.SurvEvent(cellfun(@(x) x==0, SurvEvent_RAW(:,1)),1) = {'Alive'};
 
% Convert PFI events
SURVIVAL.SurvEvent(cellfun(@(x) x==1, SurvEvent_RAW(:,2)),2) = {'Progression'};
SURVIVAL.SurvEvent(cellfun(@(x) x==0, SurvEvent_RAW(:,2)),2) = {'NoProgression'};

% Convert DSS events
SURVIVAL.SurvEvent(cellfun(@(x) x==1, SurvEvent_RAW(:,3)),3) = {'Dead'};
SURVIVAL.SurvEvent(cellfun(@(x) x==0, SurvEvent_RAW(:,3)),3) = {'Alive'};

% Convert DSS events
SURVIVAL.SurvEvent(cellfun(@(x) x==1, SurvEvent_RAW(:,4)),4) = {'Relapsed'};
SURVIVAL.SurvEvent(cellfun(@(x) x==0, SurvEvent_RAW(:,4)),4) = {'NotRelapsed'};

% Get the time columns
TimeIndx = EventIndx + 1;
SURVIVAL.SurvTime = C(:,TimeIndx);
indx_missing = ~cellfun(@(x) isnumeric(x),SURVIVAL.SurvTime);
[SURVIVAL.SurvTime{indx_missing}]  = deal(NaN);
SURVIVAL.SurvTime = cell2mat(SURVIVAL.SurvTime);


