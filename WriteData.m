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
format_str_txt = sprintf('%s%%s',Delimiter);
format_str_val = sprintf('%s%%g',Delimiter);
fprintf(fid,'Id');
if ~IdsOnly
    if  ~isempty(DATA.RowAnnotationFields)
        
        fprintf(fid,format_str_txt,DATA.RowAnnotationFields{:});
    end
end
fprintf(fid,format_str_txt,DATA.ColId{:});
% if SampleAnnotationFlag
%     fprintf(fid,'\t%s',DATA.ColId{ProbeIndx});
% end
fprintf(fid,'\n');
format_str_short = sprintf('%%s');
for i=1:DATA.nRow
    fprintf(fid,format_str_short,DATA.RowId{i});
    if ~IdsOnly
        fprintf(fid,format_str_txt,DATA.RowAnnotation{i,:});
    end
    fprintf(fid,format_str_val,DATA.X(i,:));
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
fclose(fid);






end