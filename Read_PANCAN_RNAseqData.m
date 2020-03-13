function [DATA] = Read_PANCAN_RNAseqData


%Read 'EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv'
DATA = ReadDataAEB('EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv');


% Add Info

DATA.Title = 'TCGA RNAseq PanCan';
DATA.Info.Source = 'EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv';
DATA.Info.Type = 'GeneExpression';
DATA.Info.Platform = 'RNAseqV2';


% Get Refseq Id
GeneIdSymbol = cellfun(@(x) strsplit(x,'|'),DATA.ColId,'UniformOutput',false);
GeneIdSymbol = cat(1,GeneIdSymbol{:});

DATA.ColId = GeneIdSymbol(:,2);
DATA = AddREFseqGeneInfo(DATA,[],'Replace');
DATA.ColAnnotationFields(end+1) = {'Legacy GeneSymbol'};
DATA.ColAnnotation = [DATA.ColAnnotation  GeneIdSymbol(:,1)];
 
% Add basic Sample Info based on TCHA sample Id
[SampleInfo,SampleFields] = ConvertTCGA_IdAEB(DATA.RowId);
DATA.RowAnnotation = SampleInfo;
DATA.RowAnnotationFields = SampleFields;

DATA.X(DATA.X < 0) = 0;
DATA.X = log2(DATA.X+1);


end