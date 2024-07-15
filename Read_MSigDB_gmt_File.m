function DATA = Read_MSigDB_gmt_File(InputFile)


%Read input data
[filepath,file_name,ext] = fileparts(InputFile);
[FidInputFile,message] = fopen(InputFile,'r');
if  FidInputFile == -1
    disp(InputFile)
    disp(message)
    return
end

C = textscan(FidInputFile,'%q','Delimiter','\n','TextType','string');
C = C{1};

nGeneSets = length(C);

DATA.FileName = file_name;
DATA.nGeneSets = nGeneSets;
DATA.GeneSetsNames = strings([nGeneSets 1]);
DATA.WebAddress = strings([nGeneSets 1]);
DATA.GeneIds = cell(nGeneSets,1);
DATA.nGenesInSet = zeros(nGeneSets,1);

for i=1:nGeneSets
    tmprow = C(i);
    tmp = split(tmprow);
    DATA.GeneSetsNames(i) = tmp(1);
    DATA.WebAddress(i)    = tmp(2);
    DATA.GeneIds{i}       = tmp(3:end);
    DATA.nGenesInSet(i)   = length(tmp) - 2;
end


