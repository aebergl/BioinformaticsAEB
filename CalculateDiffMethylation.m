function RESULTS_DATA = CalculateDiffMethylation(DATA,DATA_Median)


CutOff_Val = 0;
indx_nonCG = ~contains(DATA.ColId,'cg');
ChrCol=strcmp('CHR',DATA.ColAnnotationFields);
% Get X & Y probes
indx_X = strcmp('X',DATA.ColAnnotation(:,ChrCol));
indx_Y = strcmp('Y',DATA.ColAnnotation(:,ChrCol));

indx_ToUse = ~(indx_nonCG | indx_X | indx_Y);



RESULTS_DATA = CreateDataStructure(DATA.nRow,27,[],[]);

% Add Info
RESULTS_DATA.Title = 'Hyper and Hypo methylation';
RESULTS_DATA.Info.Source = inputname(1);
RESULTS_DATA.Info.CutOff_Val = CutOff_Val;

VarNames = {'Average beta value','Hyper methylation','Hypo methylation','Average Hyper methylation','Average Hypo methylation'}';

RESULTS_DATA.ColId=VarNames;
RESULTS_DATA.RowId = DATA.RowId;
RESULTS_DATA.RowAnnotationFields = DATA.RowAnnotationFields;
RESULTS_DATA.RowAnnotation = DATA.RowAnnotation;

if isfield(DATA,'SURVIVAL')
    RESULTS_DATA.SURVIVAL = DATA.SURVIVAL;
end
% 
% DATA.X=B2M(DATA.X);
% DATA_Median.X=B2M(DATA_Median.X);

median_x = median(DATA_Median.X(:,indx_ToUse),1,'omitnan');


X = bsxfun(@minus, DATA.X(:,indx_ToUse),median_x );

Chr = DATA.ColAnnotation(indx_ToUse,ChrCol);

X_pos = ones(size(X)) * NaN;
X_neg = ones(size(X)) * NaN;

indx_pos = X > CutOff_Val;
indx_neg = X < -CutOff_Val;

X_pos(indx_pos) = X(indx_pos);
X_neg(indx_neg) = X(indx_neg);


% Fin max Hyper and min Hypo

X_Lim = [ones(1,size(X,2));zeros(1,size(X,2))];
X_Lim = bsxfun(@minus, X_Lim, median_x);

X_Lim_pos = ones(size(X_Lim)) * NaN;
X_Lim_neg = ones(size(X_Lim)) * NaN;

indx_pos_max = X_Lim > CutOff_Val;
indx_neg_min = X_Lim < -CutOff_Val;


X_Lim_pos(indx_pos_max) = X_Lim(indx_pos_max);
X_Lim_neg(indx_neg_min) = X_Lim(indx_neg_min);


x_max = sum(X_Lim_pos(1,:),2,'omitnan');
x_min = sum(X_Lim_neg(2,:),2,'omitnan');


RESULTS_DATA.X(:,1) = mean(DATA.X(:,indx_ToUse),2,'omitnan');
RESULTS_DATA.X(:,2) = sum(X_pos,2,'omitnan') ./ x_max;
RESULTS_DATA.X(:,3) = sum(X_neg,2,'omitnan') ./ x_min;
RESULTS_DATA.X(:,4) = mean(X_pos,2,'omitnan');
RESULTS_DATA.X(:,5) = mean(X_neg,2,'omitnan');
ChrTxt = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'};
for i = 1:length(ChrTxt)
    indx = strcmp(ChrTxt(i),Chr);
    RESULTS_DATA.X(:,i+5) = sum(X_pos(:,indx),2,'omitnan') ./ sum(X_Lim_pos(1,indx),2,'omitnan');
end
 ChrTxt = {'Chr.01','Chr.02','Chr.03','Chr.04','Chr.05','Chr.06','Chr.07','Chr.08','Chr.09','Chr.10','Chr.11','Chr.12','Chr.13','Chr.14','Chr.15','Chr.16','Chr.17','Chr.18','Chr.19','Chr.20','Chr.21','Chr.22'}';
ChrTxt=strcat({'Hyper '},ChrTxt);
 RESULTS_DATA.ColId = [RESULTS_DATA.ColId; ChrTxt];