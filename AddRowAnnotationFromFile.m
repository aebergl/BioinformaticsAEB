function DATA = AddRowAnnotationFromFile(DATA,FileName,varargin)
% USAGE:
%   DATA = AddRowAnnotationFromFile(DATA,FileName,ColumnsToUse,Delimiter)
%   Add row annotation from file to DATA
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
        
    else
        error('%s not a valid input option',varargin{i})
    end
end

% Sample Ids to use from the DATA structure
if isempty(DATA_IdName)
    DATA_Id = DATA.RowId;
else
    indx = strcmp(DATA_IdName,DATA.RowAnnotationFields);
    if any(indx)
        DATA_Id = DATA.RowAnnotation(:,indx);
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
    opts = detectImportOptions(FileName,'Sheet',SheetName,'VariableNamingRule',VariableNamingRule,'ReadVariableNames',true);
catch
    opts = detectImportOptions(FileName,'FileType','text','Delimiter',Delimiter,'VariableNamingRule',VariableNamingRule,'ReadVariableNames',true);
end

%Select variables to import
if isempty(ColumnsToAdd)
    SelectedVariables = opts.SelectedVariableNames;
else
    [SelectedVariables]  = intersect(opts.VariableNames,ColumnsToAdd,'Stable');
    opts.SelectedVariableNames = SelectedVariables;
end

opts.VariableTypes(:) = {'char'};
% get File Id to use
if isempty(File_IdName)
    File_IdColumn = 1;
else
    File_IdColumn = find(strcmp(File_IdName, opts.SelectedVariableNames));  
end
%opts = setvartype(opts,opts.SelectedVariableNames,'char');
switch lower(fExt)
    case {'.xls','.xlsb','.xlsm','.xlsx','.xltm','.xltx'}
        if isempty(SheetName)
            C = readcell(FileName,opts,'DatetimeType','text');
        else
            C = readcell(FileName,opts,'Sheet',SheetName,'DatetimeType','text');
        end
    case {'.txt','.tsv','.csv'}
        C = readcell(FileName,opts,'DatetimeType','text');
    otherwise  % Under all circumstances SWITCH gets an OTHERWISE!
        error('Unexpected file extension: %s', fExt);
end

% Fix numeric to str
indx_numeric = cellfun(@(x) isnumeric(x),C);
C(indx_numeric) = cellfun(@(x) num2str(x),C(indx_numeric),'UniformOutput',false);

File_Id = C(:,File_IdColumn);
SelectedVariables(File_IdColumn) = [];
C(:,File_IdColumn) = [];

indx_missing = ~cellfun(@(x) ischar(x),C);
[C{indx_missing}] = deal('NA');

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
Annotation = cell(DATA.nRow,size(C,2));
Annotation(:) = {'---'};

Annotation(indx_DATA,:) = C(indx_File,:);

switch lower(AddReplace)
    case 'replace'
        DATA.RowAnnotation = Annotation;
        DATA.RowAnnotationFields = SelectedVariables ;
    case 'add'
        DATA.RowAnnotation = [DATA.RowAnnotation Annotation];
        DATA.RowAnnotationFields = [DATA.RowAnnotationFields' SelectedVariables];
end

end







