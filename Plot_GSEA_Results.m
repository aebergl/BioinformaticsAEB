function fh = Plot_GSEA_Results(File_pos,File_neg,X_value_name,Size_value_Name,CollectionType,Nmax)

NameField = "NAME";

xVal_pos = [];
SizeVal_pos = [];
YtickLabelTxt_pos = [];

xVal_neg = [];
SizeVal_neg = [];
YtickLabelTxt_neg = [];

% Read Files
%Read input data
if ~ isempty(File_pos)
    opts_pos = detectImportOptions(File_pos,'FileType','delimitedtext','VariableNamingRule','preserve',"TextType","string","Delimiter","\t")
    T_pos = readtable(File_pos,opts_pos);
    xVal_pos = table2array(T_pos(:,'X_value_name'));
    SizeVal_pos = table2array(T_pos(:,'Size_value_Name'));
    YtickLabelTxt_pos = table2array(T_pos(:,NameField));
end

if ~ isempty(File_neg)
    opts_neg=detectImportOptions(File_neg,'FileType','delimitedtext','VariableNamingRule','preserve',"TextType","string","Delimiter","\t")
    T_neg=readtable(File_neg,opts_pos);
    xVal_neg = table2array(T_neg(:,'X_value_name'));
    SizeVal_neg = table2array(T_neg(:,'Size_value_Name'));
    YtickLabelTxt_neg = table2array(T_neg(:,NameField));
end

