function RESULTS_DATA = CalculateDiffMethylation(DATA,DATA_Median)


CutOff_Val = 0;
indx_nonCG = ~contains(DATA.ColId,'cg');
ChrCol=strcmp('CHR',DATA.ColAnnotationFields);
% Get X & Y probes
indx_X = strcmp('X',DATA.ColAnnotation(:,ChrCol));
indx_Y = strcmp('Y',DATA.ColAnnotation(:,ChrCol));

indx_ToUse = ~(indx_nonCG | indx_X | indx_Y);


if isempty(DATA_Median)
    RESULTS_DATA = CreateDataStructure(DATA.nRow,23,[],[]);
    VarNames = {'Average beta value'}';

else
    RESULTS_DATA = CreateDataStructure(DATA.nRow,5+22+22,[],[]);
    VarNames = {'Average beta value','Hyper methylation','Hypo methylation','Average Hyper methylation','Average Hypo methylation'}';

end


% Add Info
RESULTS_DATA.Title = 'Hyper and Hypo methylation';
RESULTS_DATA.Info.Source = inputname(1);
RESULTS_DATA.Info.CutOff_Val = CutOff_Val;
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

if isempty(DATA_Median)
    median_x = zeros(1,sum(indx_ToUse));
else
    median_x = median(DATA_Median.X(:,indx_ToUse),1,'omitnan');
end

% Calculate X matrix with median removed
X = bsxfun(@minus, DATA.X(:,indx_ToUse),median_x );

Chr = DATA.ColAnnotation(indx_ToUse,ChrCol);

X_pos = ones(size(X)) * NaN;
X_neg = ones(size(X)) * NaN;

% Find index for hype and hypo methylayted CpG-probes
indx_pos = X > CutOff_Val;
indx_neg = X < -CutOff_Val;

% Put values into matrix
X_pos(indx_pos) = X(indx_pos);
X_neg(indx_neg) = X(indx_neg);


% Fin max Hyper and min Hypo using a samples with all ones or zeros
% these wil be used to normalise the final values

% Two extreme samples
X_Lim = [ones(1,size(X,2));zeros(1,size(X,2))]; 

% Remove median
X_Lim = bsxfun(@minus, X_Lim, median_x);

% Generate results matrix with all NaNs
X_Lim_val = ones(size(X_Lim)) * NaN;

% Find values above and below cut off
indx_pos_max = X_Lim(1,:) > CutOff_Val;
indx_neg_min = X_Lim(2,:) < -CutOff_Val;

% Put values into matrix
X_Lim_val(1,indx_pos_max) = X_Lim(1,indx_pos_max);
X_Lim_val(2,indx_neg_min) = X_Lim(2,indx_neg_min);

% Fin the maximum hyper and minimum hypo values
x_max = sum(X_Lim_val(1,:),2,'omitnan');
x_min = sum(X_Lim_val(2,:),2,'omitnan');

RESULTS_DATA.X(:,1) = mean(DATA.X(:,indx_ToUse),2,'omitnan');
ChrTxt = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'};

if isempty(DATA_Median)
    for i = 1:length(ChrTxt)
        indx = strcmp(ChrTxt(i),Chr);
        RESULTS_DATA.X(:,i+1) = mean(X(:,indx),2,'omitnan');
    end
    ChrTxt = {'Chr.01','Chr.02','Chr.03','Chr.04','Chr.05','Chr.06','Chr.07','Chr.08','Chr.09','Chr.10','Chr.11','Chr.12','Chr.13','Chr.14','Chr.15','Chr.16','Chr.17','Chr.18','Chr.19','Chr.20','Chr.21','Chr.22'}';
    ChrTxt=strcat({'Average '},ChrTxt);
else
    RESULTS_DATA.X(:,2) = sum(X_pos,2,'omitnan') ./ x_max;
    RESULTS_DATA.X(:,3) = sum(X_neg,2,'omitnan') ./ x_min;
    RESULTS_DATA.X(:,4) = mean(X_pos,2,'omitnan');
    RESULTS_DATA.X(:,5) = mean(X_neg,2,'omitnan');
    for i = 1:length(ChrTxt)
        indx = strcmp(ChrTxt(i),Chr);
        RESULTS_DATA.X(:,i+5) = sum(X_pos(:,indx),2,'omitnan') ./ sum(X_Lim_val(1,indx),2,'omitnan');
        RESULTS_DATA.X(:,i+5+22) = sum(X_neg(:,indx),2,'omitnan') ./ sum(X_Lim_val(2,indx),2,'omitnan');
    end
    ChrTxt = {'Chr.01','Chr.02','Chr.03','Chr.04','Chr.05','Chr.06','Chr.07','Chr.08','Chr.09','Chr.10','Chr.11','Chr.12','Chr.13','Chr.14','Chr.15','Chr.16','Chr.17','Chr.18','Chr.19','Chr.20','Chr.21','Chr.22'}';
    ChrTxtHyper=strcat({'Hyper '},ChrTxt);
    ChrTxtHypo=strcat({'Hypo '},ChrTxt);
end


RESULTS_DATA.ColId = [RESULTS_DATA.ColId; ChrTxtHyper; ChrTxtHypo];