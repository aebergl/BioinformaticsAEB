function [SampleInfo, ColumnNames] = Geo2SampleInfo(GEO,FileOut)

Samplefields = fieldnames(GEO.Header.Samples);
numSamplefields = length(Samplefields);
S = struct2cell(GEO.Header.Samples);
numSamples = length(S{1});

tmp = cellfun(@(x) size(x,1), S, 'UniformOutput', false);
tmp=cell2mat(tmp);
numTotalSampleFields=sum(tmp);


SampleInfo = cell(numSamples,numTotalSampleFields);
SampleInfo(:) = {'---'};
ColumnNames = cell(numTotalSampleFields,1);

counter = 0;
for i = 1:numSamplefields
    Val = S{i}';
    nVal = size(Val,2);

    if nVal == 1
        counter = counter + 1;
        ColumnNames(counter) = Samplefields(i);
        SampleInfo(:,counter) = Val(1:numSamples);
    else
        for j=1:nVal
            AddNameFlag = false;
            counter = counter + 1;
            for k=1:numSamples
                Val_1 = Val{k,j};
                if ~isempty(Val_1)
                    indx=strfind(Val_1,':');
                    if ~AddNameFlag
                        ColumnNames{counter} = Val_1(1:indx-1);
                        AddNameFlag = true;
                    end
                    Val_1 =  Val_1(indx+1:end);
                    Val_1 = strtrim(Val_1);
                    SampleInfo{k,counter} = Val_1;
                end
            end
        end
    end
end

SameValIndx = false(numTotalSampleFields,1);
for i = 1:numTotalSampleFields
    %unique(SampleInfo(:,i))
    if length(unique(SampleInfo(:,i))) == 1
        SameValIndx(i) = true;
    end
end

ColumnNames(SameValIndx) = [];
SampleInfo(:,SameValIndx) = [];

if ~isempty(FileOut)
    [fid,message] = fopen(FileOut,'w');
    if  fid == -1
        disp(FileOut)
        disp(message)
        return
    end
    fprintf(fid,'%s\t',ColumnNames{1:end-1});
    fprintf(fid,'%s\n',ColumnNames{end});
    for i=1:numSamples
        fprintf(fid,'%s\t',SampleInfo{i,1:end-1});
        fprintf(fid,'%s\n',SampleInfo{i,end});
        
    end
    fclose(fid);
end


