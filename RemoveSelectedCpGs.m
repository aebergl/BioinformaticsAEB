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

TotalRemoved = 0;

fprintf('\n')
% Only Keep cg probes
if KeepCgOnly
    indx_cg = strncmpi("cg",DATA.ColId,2);
    if any(indx_cg)
        DATA =  EditVariablesDATA(DATA,DATA.ColId(indx_cg),'Keep');
        fprintf('%u non-cg variables removed\n',sum(~indx_cg))
    end

end

% Remove chromosome X probes
if RemoveX
    indx_ChrCol = ismember(DATA.ColAnnotationFields,ChrColumnId);
    ChrData = DATA.ColAnnotation(:,indx_ChrCol);
    indx_X = contains(ChrData,ChrX);
    indx_X = any(indx_X,2);
    if any(indx_X)
        DATA =  EditVariablesDATA(DATA,DATA.ColId(indx_X),'Remove');
        fprintf('%u ChrX  variables removed\n',sum(indx_X))
        TotalRemoved = TotalRemoved + sum(indx_X);
    end
end

% Remove chromosome Y probes
if RemoveY
    indx_ChrCol = ismember(DATA.ColAnnotationFields,ChrColumnId);
    ChrData = DATA.ColAnnotation(:,indx_ChrCol);
    indx_Y = contains(ChrData,ChrY);
    indx_Y = any(indx_Y,2);
    if any(indx_Y)
        DATA =  EditVariablesDATA(DATA,DATA.ColId(indx_Y),'Remove');
        fprintf('%u ChrY variables removed\n',sum(indx_Y))
        TotalRemoved = TotalRemoved + sum(indx_Y);
    end
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
    if numel(CpG_remove) > 1
        DATA =  EditVariablesDATA(DATA,CpG_remove,'Remove');
        fprintf('%u MASKED variables removed\n',numel(CpG_remove))
        TotalRemoved = TotalRemoved + numel(CpG_remove);
    end
end

if RemoveMV
    PercMV = sum(isnan(DATA.X),1) / DATA.nRow * 100;
    indx_MV = PercMV > PercMVtol;
    if any(indx_MV)
        DATA =  EditVariablesDATA(DATA,DATA.ColId(indx_MV),'Remove');
        fprintf('%u variables removed with more than %u%% Missing Values\n',sum(indx_MV),PercMVtol)
        TotalRemoved = TotalRemoved + sum(indx_MV)
    end
end

fprintf('%u variables removed in TOTAL\n',TotalRemoved)

end


