function DATA = AddRowAnnotationFromFile(DATA,DATAIdToUse,FileName,FileIdToUse,ColumnsToAdd,SheetName)
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


[fPath, fName, fExt] = fileparts(FileName);

try 
    opts = detectImportOptions(FileName);   
catch
    opts = detectImportOptions(FileName,'FileType','text');
end

% find matching columns

indxID = strcmp(FileIdToUse,opts.VariableNames);


switch lower(fExt)
  case '.xls' || '.xlsb' || '.xlsm' || '.xlsx' || '.xltm' || '.xltx'
    C = readcell(FileName,'Sheet',SheetName,'NumHeaderLines',0);
  case '.txt'
    % A Text file
  otherwise  % Under all circumstances SWITCH gets an OTHERWISE!
    error('Unexpected file extension: %s', fExt);
end