function [LUMP] = CalcLUMP(DATA)
% High LUMP = pure sample
% Low LUMP = infiltrating immune cells


load('/Users/bergluae/AEBERGL/DATA/METHYLATION/LUMP/probes_lump.mat','probes_lump');
[~,ia,ib] = intersect(DATA.ColId,probes_lump,'stable');

% DATA.X(DATA.X == 0) = NaN;


if ~length(ia)==44 || ~length(ib)==44
    disp('Not All 44 probes found');
    
end

X = DATA.X(:,ia);


LUMP = mean(X,2,'omitnan')/0.85;

LUMP(LUMP > 1) = 1;
% DATA.X = [DATA.X;LUMP];
% DATA.P = [DATA.P;zeros(1,length(LUMP))];
% DATA.NumProbes = DATA.NumProbes + 1; 
% DATA.ProbeId(end+1) = {'LUMP'};
% DATA.ProbeAnnotation(end+1,:) = {'LUMP'};
% DATA.ProbeAnnotationColumns(end+1) = {'LUMP'};



end