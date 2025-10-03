function [DATA] = RemoveSelectedCpGs(DATA,varargin)
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

KeepCgOnly      = false;
RemoveMasked    = false;
RemoveX         = false;
RemoveY         = false;
RemoveMV        = false;
MaskedColToUse  = {'MASK_general','M_general'};
PercMVtol         = 0;

i=0;
% Check input arguments
while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'KeepCgOnly')
        KeepCgOnly = true;
    elseif strcmpi(varargin{i},'X')
        RemoveX = true;
    elseif strcmpi(varargin{i},'Y')
        RemoveY = true;
    elseif strcmpi(varargin{i},'Masked')
        RemoveMasked = true;
        i = i + 1;
        MaskedFile = varargin{i};
    elseif strcmpi(varargin{i},'MaskedColumnToUse')
        i = i + 1;
        MaskedColToUse = varargin{i};
        elseif strcmpi(varargin{i},'RemoveMV')
        RemoveMV = true;    
        i = i + 1;
        PercMVtol = varargin{i};

    end
end

% Set column Ids to check for
ChrColumnId     = ["CHR","Chromosome_36","CpG_chrm"];
ChrX            = ["X","chrX"];
ChrY            = ["Y","chrY"];

% Only Keep cg probes
if KeepCgOnly
    indx_cg = strncmpi("cg",DATA.ColId,2);
    DATA =  EditVariablesDATA(DATA,DATA.ColId(indx_cg),'Keep');
end

% Remove chromosome X probes
if RemoveX
    indx_ChrCol = ismember(DATA.ColAnnotationFields,ChrColumnId);
    ChrData = DATA.ColAnnotation(:,indx_ChrCol);
    indx_X = contains(ChrData,ChrX);
    indx_X = any(indx_X,2);
    DATA =  EditVariablesDATA(DATA,DATA.ColId(indx_X),'Remove');
end

% Remove chromosome Y probes
if RemoveY
    indx_ChrCol = ismember(DATA.ColAnnotationFields,ChrColumnId);
    ChrData = DATA.ColAnnotation(:,indx_ChrCol);
    indx_Y = contains(ChrData,ChrY);
    indx_Y = any(indx_Y,2);
    DATA =  EditVariablesDATA(DATA,DATA.ColId(indx_Y),'Remove');
end

% Remove Masked probes
if RemoveMasked
    L = readlines(MaskedFile);
    indx = L == ""; % Remove emty rows
    L(indx) = [];
    L = split(L);
    ColHead = L(1,2:end);
    MASK = L(2:end,2:end);
    CpgIds = L(2:end,1);
    indx_MaskCol = ismember(ColHead,MaskedColToUse);
    MASK = MASK(:,indx_MaskCol);
    indx_M = contains(MASK,"TRUE");
    indx_M = any(indx_M,2);
    CpG_remove = CpgIds(indx_M);
    DATA =  EditVariablesDATA(DATA,CpG_remove,'Remove');
end

if RemoveMV
    PercMV = sum(isnan(DATA.X),1) / DATA.nRow * 100;
    indx_MV = PercMV > PercMVtol;
    DATA =  EditVariablesDATA(DATA,DATA.ColId(indx_MV),'Remove');
end



end


