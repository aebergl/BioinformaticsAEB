function DATA_out = Read_GSEA_edb_File(InputFile,ResultName,GeneSetType)



% Check if File exists
if ~isfile(InputFile)
    error('Could not open %s',InputFile)
end


% Read the file
S = readstruct(InputFile,'FileType','xml');

% Define variable names
VarName = {'SIZE','ES','NES','NOM p-val','FDR q-val','FWER p-val','RANK AT MAX','RANK SCORE'};
AttributesToKeep = ["ESAttribute" "NESAttribute" "NPAttribute" "FDRAttribute" "FWERAttribute" "RANK_AT_ESAttribute" "RANK_SCORE_AT_ESAttribute"];
[~, indx_ToKeep] = ismember(AttributesToKeep,fieldnames(S.DTG));




%Get the size of the data
nRow = length(S.DTG);
nCol = length(VarName);

% Create data structure
DATA_out = CreateDataStructure(nRow,nCol,[],[]);
DATA_out.Title = S.DTG(1).RANKED_LISTAttribute;


% Convert data into a matrix
C = struct2cell(S.DTG);
DATA_out.X(:,2:end)=squeeze(cell2mat(C([indx_ToKeep],1,:)))';

% Get number of genes in each pathway
nGenes = squeeze((C(11,1,:)));
nGenes = cellfun(@(x) split(x),nGenes,'UniformOutput',false);
nGenes = cellfun(@(x) length(x),nGenes,'UniformOutput',false);
nGenes = cell2mat(nGenes);
DATA_out.X(:,1) = nGenes;

% Get the name of the gene set
Names = squeeze((C(3,1,:)));
Names = [Names{:}]';
Names = strrep(Names,'gene_sets.gmt#','');

% Assign IDs and annotation
DATA_out.ColId = string(VarName');
if isempty(ResultName)
    DATA_out.RowId = Names;
    DATA_out.RowAnnotationFields = "Gene Set Name";
    DATA_out.RowAnnotation = Names;
else
    DATA_out.RowId = strcat(ResultName," ",Names);
    DATA_out.RowAnnotationFields = ["Gene Set Name" "Result Name"];
    tmp = strings([nRow 1]);
    tmp(:) = ResultName;
    DATA_out.RowAnnotation = [Names tmp];
end

if ~isempty (GeneSetType)
switch lower(GeneSetType)
    case 'hallmark'
        GSEA_Name ={
            'HALLMARK_ADIPOGENESIS','Adipogenesis','development';
            'HALLMARK_ALLOGRAFT_REJECTION','Allograft Rejection','immune';
            'HALLMARK_ANDROGEN_RESPONSE','Androgen Response','signaling';
            'HALLMARK_ANGIOGENESIS','Angiogenesis','development';
            'HALLMARK_APICAL_JUNCTION','Apical Junction','cellular component';
            'HALLMARK_APICAL_SURFACE','Apical Surface','cellular component';
            'HALLMARK_APOPTOSIS','Apoptosis','pathway';
            'HALLMARK_BILE_ACID_METABOLISM','Bile Acid Metabolism','metabolic';
            'HALLMARK_CHOLESTEROL_HOMEOSTASIS','Cholesterol Homeostasis','metabolic';
            'HALLMARK_COAGULATION','Coagulation','immune';
            'HALLMARK_COMPLEMENT','Complement','immune';
            'HALLMARK_DNA_REPAIR','DNA Repair','DNA damage';
            'HALLMARK_E2F_TARGETS','E2F Targets','proliferation';
            'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','EMT','development';
            'HALLMARK_ESTROGEN_RESPONSE_EARLY','Estrogen Response Early','signaling';
            'HALLMARK_ESTROGEN_RESPONSE_LATE','Estrogen Response Late','signaling';
            'HALLMARK_FATTY_ACID_METABOLISM','Fatty Acid Metabolism','metabolic';
            'HALLMARK_G2M_CHECKPOINT','G2M Checkpoint','proliferation';
            'HALLMARK_GLYCOLYSIS','Glycolysis','metabolic';
            'HALLMARK_HEDGEHOG_SIGNALING','Hedgehog Signaling','signaling';
            'HALLMARK_HEME_METABOLISM','Heme Metabolism','metabolic';
            'HALLMARK_HYPOXIA','Hypoxia','pathway';
            'HALLMARK_IL2_STAT5_SIGNALING','IL2 STAT5 Signaling','signaling';
            'HALLMARK_IL6_JAK_STAT3_SIGNALING','IL6 JAK STAT3 Signaling','immune';
            'HALLMARK_INFLAMMATORY_RESPONSE','Inflammatory Response','immune';
            'HALLMARK_INTERFERON_ALPHA_RESPONSE','Interferon alpha Response','immune';
            'HALLMARK_INTERFERON_GAMMA_RESPONSE','Interferon gamma Response','immune';
            'HALLMARK_KRAS_SIGNALING_DN','KRAS Signaling Dn','signaling';
            'HALLMARK_KRAS_SIGNALING_UP','KRAS Signaling Up','signaling';
            'HALLMARK_MITOTIC_SPINDLE','Mitotic Spindle','proliferation';
            'HALLMARK_MTORC1_SIGNALING','MTORC1 Signaling','signaling';
            'HALLMARK_MYC_TARGETS_V1','MYC Targets V1','proliferation';
            'HALLMARK_MYC_TARGETS_V2','MYC Targets V2','proliferation';
            'HALLMARK_MYOGENESIS','Myogenesis','development';
            'HALLMARK_NOTCH_SIGNALING','NOTCH Signaling','signaling';
            'HALLMARK_OXIDATIVE_PHOSPHORYLATION','Oxidative Phosphorylation','metabolic';
            'HALLMARK_P53_PATHWAY','P53 Pathway','proliferation';
            'HALLMARK_PANCREAS_BETA_CELLS','Pancreas beta Cell','development';
            'HALLMARK_PEROXISOME','Peroxisome','cellular component';
            'HALLMARK_PI3K_AKT_MTOR_SIGNALING','PI3K AKT MTOR Signaling','signaling';
            'HALLMARK_PROTEIN_SECRETION','Protein Secretion','pathway';
            'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY','Reactive Oxygen Species','pathway';
            'HALLMARK_SPERMATOGENESIS','Spermatogenesis','development';
            'HALLMARK_TGF_BETA_SIGNALING','TGF beta Signaling','signaling';
            'HALLMARK_TNFA_SIGNALING_VIA_NFKB','TNFa Signaling Via NFkB','signaling';
            'HALLMARK_UNFOLDED_PROTEIN_RESPONSE','Unfolded Protein Response','pathway';
            'HALLMARK_UV_RESPONSE_DN','UV Response Down','DNA damage';
            'HALLMARK_UV_RESPONSE_UP','UV Response Up','DNA damage';
            'HALLMARK_WNT_BETA_CATENIN_SIGNALING','WNT beta Catenin Signaling','signaling';
            'HALLMARK_XENOBIOTIC_METABOLISM','Xenobiotic Metabolism','metabolic';
            };
        [~,indx_new]=ismember(Names,GSEA_Name(:,1));

        ShortName = GSEA_Name(indx_new,2);
        Group = GSEA_Name(indx_new,3);
        DATA_out.RowAnnotationFields = horzcat(DATA_out.RowAnnotationFields,["Short Gene Set Name" "Group"]);

        DATA_out.RowAnnotation = [DATA_out.RowAnnotation ShortName Group];


    case 'reactome'
        YtickLabelTxt = strrep(YtickLabelTxt,'REACTOME_','');
        YtickLabelTxt = strrep(YtickLabelTxt,'_',' ');
        YtickLabelTxt = lower(YtickLabelTxt);
        YtickLabelTxt = upper(extractBefore(YtickLabelTxt,2)) + extractAfter(YtickLabelTxt,1);
end


end