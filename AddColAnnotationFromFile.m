function DATA = AddColAnnotationFromFile(DATA,FileName,varargin)
% USAGE:
%   DATA = AddColAnnotationFromFile(DATA,FileName,varargin)
%   Add column annotation from file to DATA
%
% INPUTS:
% * DATA: DATA structure
% * IdToUse: Row identifier to use, [], RowId is used
% * SURVIVAL: SURVIVAL structure
%
% OUTPUTS:
% * DATA: DATA structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Anders Berglund, 2020 aebergl at gmail.com                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SheetName = '';
DATA_IdName = [];
File_IdName = [];
ColumnsToAdd = [];
Truncate = 0;
AddReplace = 'Add';
Delimiter = {'\t'};
VariableNamingRule = 'preserve';
% Check Input
i=0;
while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'DATA_Id')
        i = i + 1;
        DATA_IdName = varargin{i};
    elseif strcmpi(varargin{i},'File_Id')
        i = i + 1;
        File_IdName = varargin{i};
    elseif strcmpi(varargin{i},'ColumnsToAdd')
        i = i + 1;
        ColumnsToAdd = varargin{i};
    elseif strcmpi(varargin{i},'Delimiter')
        i = i + 1;
        Delimiter = varargin{i};

    elseif strcmpi(varargin{i},'SheetName')
        i = i + 1;
        SheetName = varargin{i};
    elseif strcmpi(varargin{i},'Truncate')
        i = i + 1;
        Truncate = varargin{i};        
    elseif strcmpi(varargin{i},'Replace')
        AddReplace = 'Replace';
    elseif strcmpi(varargin{i},'Add')
        AddReplace = 'Add';        
    end
end

% Sample Ids to use from the DATA structure
if isempty(DATA_IdName)
    DATA_Id = DATA.ColId;
else
    indx = strcmp(DATA_IdName,DATA.ColAnnotationFields);
    if any(indx)
        DATA_Id = DATA.ColAnnotation(:,indx);
    else
        error('Could not find the given Id in the DATA strcuture')
    end   
end

% Get info for the annotation file
try
    [fPath, fName, fExt] = fileparts(FileName);
catch
    error('Could not find file: %s',FileName);
end

% Get info for annotation file
try
    opts = detectImportOptions(FileName,'Sheet',SheetName,'VariableNamingRule',VariableNamingRule);
catch
    opts = detectImportOptions(FileName,'FileType','text','Delimiter',Delimiter,'VariableNamingRule',VariableNamingRule);
end

%Select variables to import
% if isempty(ColumnsToAdd)
%     SelectedVariables = opts.SelectedVariableNames;
% else
%     [SelectedVariables]  = intersect(ColumnsToAdd,opts.VariableNames,'Stable');
%     opts.SelectedVariableNames = SelectedVariables;
% end

% get File Id to use
if isempty(File_IdName)
    File_IdColumn = 1;
else
    File_IdColumn = find(strcmp(File_IdName, opts.VariableNames));  
end
opts = setvartype(opts,opts.VariableNames,'char');
switch lower(fExt)
    case {'.xls','.xlsb','.xlsm','.xlsx','.xltm','.xltx'}
        if isempty(SheetName)
            C = readcell(FileName,opts);
        else
            C = readcell(FileName,opts,'Sheet',SheetName);
        end
    case {'.txt','.tsv','.csv','.annot'}
        C = readcell(FileName,opts,'TextType','char');
    otherwise  % Under all circumstances SWITCH gets an OTHERWISE!
        error('Unexpected file extension: %s', fExt);
end
 
% Fix numeric to str
indx_numeric = cellfun(@(x) isnumeric(x),C);
C(indx_numeric) = cellfun(@(x) num2str(x),C(indx_numeric),'UniformOutput',false);

% Get File Id column
File_Id = C(:,File_IdColumn);

%SelectedVariables(File_IdColumn) = [];
%C(:,File_IdColumn) = [];


% Extract slected Columns to add
[SelectedVariables,~,ib]  = intersect(ColumnsToAdd,opts.VariableNames,'Stable');
C = C(:,ib);

indx_missing = ~cellfun(@(x) ischar(x),C);
[C{indx_missing}] = deal('NA');

%Remove rows with missing ids for merging
indx_missing = ~cellfun(@(x) ischar(x),File_Id);
C(indx_missing,:) = [];
File_Id(indx_missing,:) = [];


if Truncate    
    File_Id = cellfun(@(x) x(1:Truncate), File_Id, 'UniformOutput', false);
    DATA_Id = cellfun(@(x) x(1:Truncate), DATA_Id, 'UniformOutput', false);    
end

[indx_DATA,indx_File]  = ismember(DATA_Id,File_Id);
indx_File = indx_File(indx_File>0);

%Create Annotation object
Annotation = cell(DATA.nCol,size(C,2));
Annotation(:) = {'---'};

Annotation(indx_DATA,:) = C(indx_File,:);

% Convert everyhting to string arrays
Annotation = string(Annotation);
SelectedVariables = string(SelectedVariables);

switch lower(AddReplace)
    case 'replace'
        DATA.ColAnnotation = Annotation;
        DATA.ColAnnotationFields = SelectedVariables ;
    case 'add'
        if iscellstr(DATA.ColAnnotation)
            DATA.ColAnnotation = string(DATA.ColAnnotation);
        end
       if iscellstr(DATA.ColAnnotationFields)
            DATA.ColAnnotationFields = string(DATA.ColAnnotationFields);
       end
        DATA.ColAnnotation = [DATA.ColAnnotation Annotation];
        DATA.ColAnnotationFields = [DATA.ColAnnotationFields'; SelectedVariables'];
end

end







