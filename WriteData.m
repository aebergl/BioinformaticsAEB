function WriteData(DATA,FileOut,varargin)

Delimiter = '\t';
SampleAnnotationFlag = false;
IdsOnly = false;
RowIds = [];
ColIds = [];
i=0;
while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'Delimiter')
        i = i + 1;
        Delimiter = varargin{i};
    elseif strcmpi(varargin{i},'SampleAnnotation')
        SampleAnnotationFlag = true;
    elseif strcmpi(varargin{i},'IdsOnly')
        IdsOnly = true;
    elseif strcmpi(varargin{i},'ProbeId')
        i = i + 1;
        ProbeIds = varargin{i};
        [~,~,ProbeIndx]=intersect(ProbeIds,DATA.ProbeId,'stable');
    elseif strcmpi(varargin{i},'SampleId')
        i = i + 1;
        SampleIds = varargin{i};
    end
end

[fid,message] = fopen(FileOut,'w');
if  fid == -1
    disp(FileOut)
    disp(message)
    return
end

fprintf(fid,'Id');
if  ~isempty(DATA.RowAnnotationFields)
    format_str = sprintf('%s%%s',Delimiter);
    fprintf(fid,format_str,DATA.RowAnnotationFields{:});
end
% if SampleAnnotationFlag
%     fprintf(fid,'\t%s',DATA.ColId{ProbeIndx});
% end
fprintf(fid,'\n');
format_str_short = sprintf('%%s');
for i=1:DATA.nRow
    fprintf(fid,format_str_short,DATA.RowId{i});
    fprintf(fid,format_str,DATA.RowAnnotation{i,:});
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
fclose(fid);






end