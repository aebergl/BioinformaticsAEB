function NumLines = GetNumLines(FileName)

[fid,message] = fopen(FileName,'r');
if  fid == -1
    disp(FileName)
    disp(message)
    return
end

fseek(fid, 0, 'eof');
chunksize = ftell(fid);
fseek(fid, 0, 'bof');
ch = fread(fid, chunksize, '*uchar');
NumLines = sum(ch == newline);
if NumLines == 0 
    NumLines = sum(ch == sprintf('\r')) + 1;
end
fclose(fid);
