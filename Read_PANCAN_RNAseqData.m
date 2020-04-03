function [DATA] = Read_PANCAN_RNAseqData


%Read already downloaded 'EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv'
if ~isfile('EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv')
    try
        fprintf('Could not find EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv in this directory\n')
        fprintf('Starting to download it now, will probably take a couple of minutes\n')
        outfilename = websave('EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv','http://api.gdc.cancer.gov/data/3586c0da-64d0-4b74-a449-5ff4d9136611');
        
    catch
        error('Could not load EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv from https://gdc.cancer.gov/about-data/publications/pancanatlas')
    end
end
fprintf('Processing EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv\n')

DATA = ReadData('EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv');

% Add Info
DATA.Title = 'TCGA RNAseq PanCan';
DATA.Info.Source = 'EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv';
DATA.Info.Type = 'GeneExpression';
DATA.Info.Platform = 'RNAseqV2';


% Get Refseq Id
fprintf('Updating Gene Annotation\n')

GeneIdSymbol = cellfun(@(x) strsplit(x,'|'),DATA.ColId,'UniformOutput',false);
GeneIdSymbol = cat(1,GeneIdSymbol{:});
DATA.ColId = GeneIdSymbol(:,2);
DATA = AddREFseqGeneInfo(DATA,[],'Replace');
DATA.ColAnnotationFields(end+1:end+2) = {'TCGA GeneSymbol','TCGA GeneId'}';
DATA.ColAnnotation = [DATA.ColAnnotation  GeneIdSymbol];

% Removing samples based on 'Merged Sample Quality Annotations - merged_sample_quality_annotations.tsv'
IdsToRemove = CheckTCGASampleQuality(DATA.RowId,[]);
DATA  = EditSamplesDATA(DATA,IdsToRemove,'Remove');



% Add basic Sample info based on TCGA sample Id
fprintf('Converting TCGA Id to sample info\n')

[SampleInfo,SampleFields] = ConvertTCGAId(DATA.RowId);
DATA.RowAnnotation = SampleInfo;
DATA.RowAnnotationFields = SampleFields;

% Convert to log2
DATA.X(DATA.X < 0) = 0;
DATA.X = log2(DATA.X+1);

% Add Survival data from TCGA-Clinical Data Resource (CDR) Outcome
% TCGA-CDR-SupplementalTableS1.xlsx 
% Will add Survival dat for the normal samples
fprintf('Adding Survival data\n')
SURVIVAL = ReadTCGA_Survival_File;

DATA = AddSurvivalDATA(DATA,[],SURVIVAL);


end