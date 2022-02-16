function DATA = ReadMultipleTxtFiles(ListOfFiles)
% USAGE:
%   DATA = AddSurvivalDATA(DATA,FileName,ColumnsToUse,Delimiter)
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

nFiles = length(ListOfFiles);

% Get info from the first file
numRows = GetNumLines(FileName);
opts = detectImportOptions(FileName,'FileType','text');

nRowsData = numRows - (opts.DataLines(1) - 1);

DATA = CreateDataStructure(nRowsData,nFiles,1,1);

for i=1:nFiles
    
end