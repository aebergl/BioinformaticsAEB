function NumColumn = GetNumColumns(FileName,DelimiterType,NumnLinesToRead)


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