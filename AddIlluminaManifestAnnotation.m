function DATA  = AddIlluminaManifestAnnotation(DATA,AnnotationFile,varargin)
% USAGE:
%   DATA = AddColAnnotationFromFile(DATA,FileName,varargin)
%   Add column annotation from Illumina manifest file, swork for M450k EPICv1 and EPICv2
%
% INPUTS:
% * DATA: DATA structure, CpG-probeId in DATA.ColId
% * AnnotationFile: Illumina manifest file to use, must match chip type in DATA
%
% OUTPUTS:
% * DATA: DATA structure with columns annotated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Anders Berglund, 2025 aebergl at gmail.com                           %

AddReplace      = "Add";
CpG_Id          = "IlmnID";
ColumnsToUse    = [];

i=0;
while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'Replace')
        AddReplace = 'Replace';
    elseif strcmpi(varargin{i},'Add')
        AddReplace = 'Add';
    elseif strcmpi(varargin{i},'ColumnsToAdd')
        i = i + 1;
        ColumnsToUse = varargin{i};
    end

end

% read files into a string array
L = readlines(AnnotationFile);

% Find different sections
AssayPos     = find(contains(L,'[Assay]'));
ControlsPos  = find(contains(L,'[Controls]'));
LociCountPos = contains(L,'Loci Count');

AnnotationColumnIds = split(L(AssayPos+1),',');


% Find number of CpG probes in the Annotation file
tmp = L(LociCountPos);
tmp = strsplit(tmp,',');
nLoci = str2double(tmp(2));
if nLoci < DATA.nCol
    warning('Data strucure have more columns (n=%i) than entries in the given annotation file (n=%i)',DATA.nCol,nLoci)
end

if isempty(ColumnsToUse)
    if nLoci == 485553 % M450K
        ColumnsToUse = {'IlmnID','Infinium_Design_Type','Color_Channel','Genome_Build','CHR','MAPINFO',...
            'Strand','Probe_SNPs','Probe_SNPs_10','Random_Loci','UCSC_RefGene_Name','UCSC_RefGene_Accession',...
            'UCSC_RefGene_Group','UCSC_CpG_Islands_Name','Relation_to_UCSC_CpG_Island','Phantom','DMR',...
            'Enhancer','HMM_Island','Regulatory_Feature_Name','Regulatory_Feature_Group','DHS'};
    elseif nLoci == 865918 % EPICv1
        ColumnsToUse = {'IlmnID','Infinium_Design_Type','Genome_Build','CHR','MAPINFO',...
            'UCSC_RefGene_Name','UCSC_RefGene_Accession','UCSC_RefGene_Group','UCSC_CpG_Islands_Name',...
        	'Relation_to_UCSC_CpG_Island','450k_Enhancer','HMM_Island','Regulatory_Feature_Name','Regulatory_Feature_Group',...
            'GencodeBasicV12_NAME','GencodeBasicV12_Accession','GencodeBasicV12_Group','GencodeCompV12_NAME','GencodeCompV12_Accession','GencodeCompV12_Group',...
            'DNase_Hypersensitivity_Evidence_Count','OpenChromatin_NAME','OpenChromatin_Evidence_Count','TFBS_NAME','TFBS_Evidence_Count',...
            'SNP_ID','SNP_DISTANCE','SNP_MinorAlleleFrequency',...
            'CHR_hg38','Start_hg38','End_hg38','Strand_hg38'};
    elseif nLoci == 937055 % EPICv2
        ColumnsToUse = {'IlmnID','Name','Infinium_Design_Type','CHR','MAPINFO','Genome_Build',...
            'UCSC_RefGene_Group','UCSC_RefGene_Name','UCSC_RefGene_Accession',...
            'UCSC_CpG_Islands_Name','Relation_to_UCSC_CpG_Island',...
            'GencodeV41_Group','GencodeV41_Name','GencodeV41_Accession'...
            'Phantom5_Enhancers','HMM_Island','Regulatory_Feature_Name','Regulatory_Feature_Group'...
        	'450k_Enhancer','DMR','DNase_Hypersensitivity_NAME','Encode_CisReg_Site','Encode_CisReg_Site_Evid'...
        	'OpenChromatin_NAME','OpenChromatin_Evidence_Count',...
            'Methyl450_Loci','Methyl27_Loci','EPICv1_Loci'...
        	'Manifest_probe_match','SNP_ID','SNP_DISTANCE','SNP_MinorAlleleFrequency'};
    end
end

%Select Columns to use
[~,indx_ColumnsToUse,indx_Annotation] = intersect(ColumnsToUse,AnnotationColumnIds,'stable');
if length(indx_ColumnsToUse) < length(ColumnsToUse)
    warning('Not all columns to use fund in given annotation file')
    ColumnsToUse = ColumnsToUse(indx_ColumnsToUse);
end

% Create Annotation string array
Annotation = strings(DATA.nCol,length(indx_ColumnsToUse));
Annotation(:) = "---";

% Select block with annotation
L = L(AssayPos+2:ControlsPos-1);

%Split and create string array matrix
L = split(L,',');

% Get Cpg Id to use for merging
Annotation_CpgId = L(:,contains(AnnotationColumnIds,CpG_Id));

% Select Columns to use
L = L(:,indx_Annotation);

%Select matching CpG probes in DATA file
[indx_DATA,indx_Annotation]  = ismember(DATA.ColId,Annotation_CpgId);
indx_Annotation = indx_Annotation(indx_Annotation>0);
Annotation(indx_DATA,:) = L(indx_Annotation,:);

if any(~indx_DATA)
    warning('Not all ColIds fund in given annotation file. Some Columns will be missing annotation','backtrace','off')
end

if isempty(DATA.ColAnnotation)
    AddReplace = 'replace';
end

switch lower(AddReplace)
    case 'add'
        if iscellstr(DATA.ColAnnotation)
            DATA.ColAnnotation = string(DATA.ColAnnotation);
        end
        if iscellstr(DATA.ColAnnotationFields)
            DATA.ColAnnotationFields = string(DATA.ColAnnotationFields);
        end
        DATA.ColAnnotation = [DATA.ColAnnotation Annotation];
        DATA.ColAnnotationFields(end+1:end+length(ColumnsToUse)) = ColumnsToUse;
    case 'replace'
        DATA.ColAnnotation = Annotation;
        DATA.ColAnnotationFields = ColumnsToUse ;
end

