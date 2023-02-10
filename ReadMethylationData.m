function DATA = ReadMethylationData(InputFile,varargin)
%DATA = ReadData('GSE19188_tumor_iron_log2.txt','HeaderRows',1,'IdColumns',1)
%  reads in a text file and creates a microarray data structure
%  First row and column must be a probe/sample id
%
% DATA = ReadDataAEB(InputFile,varargin)
%
%   INPUT
%       InputFile       File with Gene expression data

%   options ---------------------------------------
%
%   'HeaderRows'            Number of Header Rows to read (default 1)
%   'HeaderRowsToIgnore'    Number of Header Rows to be ignored (default 0)
%   'EndRows'               Number of End Rows to be ignored (default 0)
%   'IdColumns'             Number of Starting Id columns (default 1)
%   'BetaValueExt'          Extension used for beta-value (default '_b')
%   'DetectionPvalExt'      Extension used for detection p-value (default '_p')
%   'T'                     Transposes the input file [DEFAULT]
%   'NoT'                   Do not transpose the input file
%   'Delimiter'             Delimiter type (default '\t')
%   'R'                     First element is missing in top row
%
% Anders Berglund
numHeaderRowsToIgnore = 0;
numHeaderRows = 1;
numEndRows = 0;
numIdColumns = 1;
DelimiterType = '\t';
TransposeX = true;
BetaValueExt = '_b';
DetectionPvalExt = '_p';
DATA=[];
R_Input = true;
Detection_Pval_CutOff = 0.01;

% Check Input
i=0;
while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'HeaderRows')
        i = i + 1;
        numHeaderRows = varargin{i};
    elseif strcmpi(varargin{i},'EndRows')
        i = i + 1;
        numEndRows = varargin{i};
    elseif strcmpi(varargin{i},'IdColumns')
        i = i + 1;
        numIdColumns = varargin{i};
    elseif strcmpi(varargin{i},'BetaValueExt')
        i = i + 1;
        BetaValueExt = varargin{i};
    elseif strcmpi(varargin{i},'DetectionPvalExt')
        i = i + 1;
        DetectionPvalExt = varargin{i};
    elseif strcmpi(varargin{i},'Delimiter')
        i = i + 1;
        DelimiterType = varargin{i};      
    elseif strcmpi(varargin{i},'T')        
        TransposeX = true;
   elseif strcmpi(varargin{i},'NoT')
        TransposeX = false;
    elseif strcmpi(varargin{i},'R')
        R_Input = 1;
    end
end
% Get info about File
numRows = GetNumLines(InputFile);
numColumns = GetNumColumns(InputFile,DelimiterType,min([numRows, 10]));


if numHeaderRows
    HeaderData = cell(numHeaderRows,1);
end

numDataRows = numRows - numHeaderRows - numEndRows - numHeaderRowsToIgnore;
numXColumns = (numColumns-numIdColumns);

%Read input data
[~,file_name,~] = fileparts(InputFile);
[FidInputFile,message] = fopen(InputFile,'r');
if  FidInputFile == -1
    warning(InputFile)
    warning(message)
    return
end

% Read top lines to be ignored
for i = 1:numHeaderRowsToIgnore
    [~] = fgetl(FidInputFile);
end

%Read Header data
for i = 1:numHeaderRows
    HeaderData{i} = fgetl(FidInputFile);
end


%Read Methylation data
ReadFormat = strcat(repmat('%q',1,numIdColumns),repmat('%f',1,numXColumns));
[S] = textscan(FidInputFile,ReadFormat,'delimiter',DelimiterType,'TreatAsEmpty',{'NA','na'});

fclose(FidInputFile);

numRead =numel(S{2});

if numRead < numDataRows

    warning('WARNING!!! %u lines read but there should have been %u',numRead,numDataRows);
    
end

DATA = CreateDataStructure(numRead,numXColumns/2,numIdColumns,numHeaderRows);


if numHeaderRows > 0
    [tmp,~] = textscan(HeaderData{1},'%q','delimiter',DelimiterType);
    tmp = tmp{1};
    if R_Input
        tmp = [{'SampleId'};tmp];
    end
    beta_indx = strfind(tmp,BetaValueExt);
    BetaAvg_Columns = (~cellfun('isempty',beta_indx));
    pval_indx = strfind(tmp,DetectionPvalExt);
    PVal_Columns = (~cellfun('isempty',pval_indx));

    tmp_names=tmp(BetaAvg_Columns);
    for i=1:numel(tmp_names)
        tmp_txt=strsplit(tmp_names{i},BetaValueExt);
        tmp_txt =   deblank(tmp_txt{1});
        DATA.ColId(i)=tmp_txt;
    end
    %DATA.ColId = string(DATA.ColId);
    if numIdColumns > 1
        DATA.RowAnnotationFields = tmp(2:numIdColumns);
    end
    for i=2:numHeaderRows
        [tmp,~] = textscan(HeaderData{i},'%q','delimiter',DelimiterType);
        tmp = tmp{1};
        DATA.ColAnnotationFields(i) = tmp(1,1);
        DATA.ColAnnotation(:,i) = tmp(numIdColumns+1:end,1);
    end
    DATA.RowId  = string(S{1});
    if numIdColumns > 1
        DATA.RowAnnotation  = [S{2:numIdColumns}];

    end
else

    
end
P = cell2mat(S(PVal_Columns));

DATA.X = cell2mat(S(BetaAvg_Columns));


DATA.X(P >= Detection_Pval_CutOff) = NaN;
DATA.X(isempty(DATA.X)) = NaN;

DATA.Title = file_name;
DATA.ShortTitle = file_name;

DATA.Info.Source = InputFile;
DATA.Info.Type = 'Methylation';


    

% Check for unique identifiers

% Start with Row Ids
numUniqueRowIds = length(unique(DATA.RowId));
if numUniqueRowIds < length(DATA.RowId)
    fprintf('WARNING!!! Row Ids are not unique\n');
end
% and now to Column Ids
numUniqueColIds = length(unique(DATA.ColId));
if numUniqueColIds < length(DATA.ColId)
    fprintf('WARNING!!! Column Ids are not unique\n');
end

% transpose DATA
if TransposeX
    DATA = TransposeData(DATA);  
end



