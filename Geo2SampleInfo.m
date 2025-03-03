function [SampleInfo, FieldNames, HeaderInfo] = Geo2SampleInfo(GEO,FileOut)
% USAGE:
%   SampleInfo, ColumnNames, HeaderInfo] = Geo2SampleInfo(GEO,FileOut)
%   process a GEO file ro readable data
%
% INPUTS:
%   DATA:       GEO structure
%   FileOut:    FilenName for output results, [] will not write the results to a file
% OUTPUTS:
%   SampleInfo: String array with sample info
%   FieldNames: String arrays with extrated fields
%   HeaderInfo: Information strucure about the series
%
%   options ---------------------------------------
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Anders Berglund, 2025 aebergl at gmail.com                           %

% Series info
HeaderInfo = GEO.Header.Series;

Samplefields = fieldnames(GEO.Header.Samples);
numSamplefields = length(Samplefields);
S = struct2cell(GEO.Header.Samples);
numSamples = length(S{1});

tmp = cellfun(@(x) size(x,1), S, 'UniformOutput', false);
tmp=cell2mat(tmp);
numTotalSampleFields=sum(tmp);


SampleInfo = cell(numSamples,numTotalSampleFields);
SampleInfo(:) = {'---'};
FieldNames = cell(numTotalSampleFields,1);

counter = 0;
for i = 1:numSamplefields
    Val = S{i}';
    nVal = size(Val,2);

    if nVal == 1
        counter = counter + 1;
        FieldNames(counter) = Samplefields(i);
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
                        FieldNames{counter} = Val_1(1:indx-1);
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

%Remove all fields that have the same value
SameValIndx = false(numTotalSampleFields,1);
for i = 1:numTotalSampleFields
    if length(unique(SampleInfo(:,i))) == 1
        SameValIndx(i) = true;
    end
end

FieldNames(SameValIndx) = [];
SampleInfo(:,SameValIndx) = [];

%Write results to a text file
if ~isempty(FileOut)
    [fid,message] = fopen(FileOut,'w');
    if  fid == -1
        disp(FileOut)
        disp(message)
        return
    end
    fprintf(fid,'%s\t',FieldNames{1:end-1});
    fprintf(fid,'%s\n',FieldNames{end});
    for i=1:numSamples
        fprintf(fid,'%s\t',SampleInfo{i,1:end-1});
        fprintf(fid,'%s\n',SampleInfo{i,end});
        
    end
    fclose(fid);
end


