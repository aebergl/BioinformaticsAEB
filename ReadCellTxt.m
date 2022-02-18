function S = ReadCellTxt(FileName)

Delimiter = '\t';

[FidInputFile,message] = fopen(FileName,'r');
if  FidInputFile == -1
    disp(RefseqGeneInfoFile)
    disp(message)
    return
end

% Get info about File
numRows = GetNumLines(FileName);
numColumns = GetNumColumns(FileName,Delimiter,min([numRows, 10]));

S = textscan(FidInputFile,repmat('%s',1,numColumns),'delimiter',Delimiter);

S = cat(2,S{1:end});

