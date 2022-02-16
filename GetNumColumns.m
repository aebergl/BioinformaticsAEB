function NumColumn = GetNumColumns(FileName,DelimiterType,NumnLinesToRead)
% NumColumn = GetNumColumns(FileName,DelimiterType,NumnLinesToRead)
%   Get number of coumns in a file
% Input
%   FileName:           Name of file
%   DelimiterType:      Delimiter to use
%   NumnLinesToRead:    Number of lines to use
% Output
%   NumColumns:         Maximum number of column found

[fid,message] = fopen(FileName,'r');
if  fid == -1
    disp(FileName)
    disp(message)
    return
end
NumColumnRow = zeros(NumnLinesToRead,1);
for i = 1:NumnLinesToRead
    tline = fgetl(fid);
    if ~isempty(tline)
        tmp = textscan(tline,'%s','delimiter',DelimiterType);
        tmp=tmp{1};
        NumColumnRow(i) = numel(tmp);
    end
end
NumColumn = max(NumColumnRow);