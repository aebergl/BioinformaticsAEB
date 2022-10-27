function RESULTS_DATA = CalculateDiffMethylation(DATA,DATA_Median)


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

VarNames = {'Average beta value','Hyper methylation','Hypo methylation','Average Hyper methylation','Average Hypo methylation'}';

RESULTS_DATA.ColId=VarNames;
RESULTS_DATA.RowAnnotationFields = DATA.RowAnnotationFields;
RESULTS_DATA.RowAnnotation = DATA.RowAnnotation;

if isfield(DATA,'SURVIVAL')
    RESULTS_DATA.SURVIVAL = DATA.SURVIVAL;

end

DATA.X=B2M(DATA.X);
DATA_Median.X=B2M(DATA_Median.X);


X = bsxfun(@minus, DATA.X(:,indx_ToUse), median(DATA_Median.X(:,indx_ToUse),1,'omitnan'));

Chr = DATA.ColAnnotation(indx_ToUse,ChrCol);

X_pos = ones(size(X)) * NaN;
X_neg = ones(size(X)) * NaN;

indx_pos = X > 0;
indx_neg = X < 0;

X_pos(indx_pos) = X(indx_pos);
X_neg(indx_neg) = X(indx_neg);

RESULTS_DATA.X(:,1) = mean(DATA.X(:,indx_ToUse),2,'omitnan');
RESULTS_DATA.X(:,2) = sum(X_pos,2,'omitnan');
RESULTS_DATA.X(:,3) = sum(X_neg,2,'omitnan');
RESULTS_DATA.X(:,4) = RESULTS_DATA.X(:,2) ./ sum(indx_pos,2);
RESULTS_DATA.X(:,5) = RESULTS_DATA.X(:,3) ./ sum(indx_neg,2);
ChrTxt = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'};
for i = 1:length(ChrTxt)
    indx = strcmp(ChrTxt(i),Chr);
    RESULTS_DATA.X(:,i+5) = sum(X_pos(:,indx),2,'omitnan');
end

