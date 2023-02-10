function [PosVal, NegVal] = CalcVolcanoStatResults(DATA,X_Variable,Y_Variable,varargin)

% Select X value
indx_X_Val = strcmpi(X_Variable,DATA.ColId);
if any(indx_X_Val)
    x_data = DATA.X(:,indx_X_Val);
else
    error('%s not found',X_Variable)
end
switch X_Variable
    case 'Delta Average'
       % x_data = x_data;
    case 'HR logrank DSS'
        x_data =  log2(x_data);
    case 'HR coxreg DSS'
        x_data =  log2(x_data);
    case 'HR logrank PFI'
        x_data =  log2(x_data);
    case 'HR coxreg PFI'
        x_data =  log2(x_data);
    case 'HR logrank OS'
        x_data =  log2(x_data);
    case 'HR coxreg OS'
        x_data =  log2(x_data);

end

% Select Y value
indx_Y_Val = strcmpi(Y_Variable,DATA.ColId);
if any(indx_Y_Val)
    y_data = DATA.X(:,indx_Y_Val);
else
    error('%s not found',Y_Variable)
end
switch Y_Variable
    case 'p logrank PFI'
        y_data =  -log10(y_data);
    case 'q logrank PFI'
        y_data =  -log10(y_data);
    case 'fdr logrank PFI'
        y_data =  -log10(y_data);
    case 'p coxreg PFI'
        y_data =  -log10(y_data);
    case 'q coxreg PFI'
        y_data =  -log10(y_data);
    case 'fdr coxreg PFI'
        y_data =  -log10(y_data);
    case 'p logrank DSS'
        y_data =  -log10(y_data);
    case 'q logrank DSS'
        y_data =  -log10(y_data);
    case 'fdr logrank DSS'
        y_data =  -log10(y_data);
    case 'p coxreg DSS'
        y_data =  -log10(y_data);
    case 'q coxreg DSS'
        y_data =  -log10(y_data);
    case 'fdr coxreg DSS'
        y_data =  -log10(y_data);
    case 'p logrank OS'
        y_data =  -log10(y_data);
    case 'q logrank OS'
        y_data =  -log10(y_data);
    case 'fdr logrank OS'
        y_data =  -log10(y_data);
    case 'p coxreg OS'
        y_data =  -log10(y_data);
    case 'q coxreg OS'
        y_data =  -log10(y_data);
    case 'fdr coxreg OS'
        y_data =  -log10(y_data);
    case 'p t-test'
        y_data =  -log10(y_data);
end

indx_pos = x_data > 0;
indx_neg = x_data < 0;

x_data = normalize(x_data,'scale');
y_data = normalize(y_data,'scale');



TotalSum = sum(sqrt(x_data.^2 + y_data.^2),'omitnan');
PosSum = sum(sqrt(x_data(indx_pos).^2 + y_data(indx_pos).^2),'omitnan');
NegSum = sum(sqrt(x_data(indx_neg).^2 + y_data(indx_neg).^2),'omitnan');

PosVal = (PosSum-NegSum)/TotalSum;
NegVal = NegSum/TotalSum;


