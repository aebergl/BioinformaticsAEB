function WriteData(DATA,FileOut,varargin)

Delimiter = '\t';
DATAFlag = true;
SampleAnnotationFlag = true;
SurvivalFlag = false;
IdsOnly = false;
SeperateFiles = false;
RowIds = [];
ColIds = [];
i=0;
while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'Delimiter')
        i = i + 1;
        Delimiter = varargin{i};
    elseif strcmpi(varargin{i},'SampleAnnotationOnly')
        SampleAnnotationFlag = true;
        DATAFlag = false;
    elseif strcmpi(varargin{i},'Survival')
        SurvivalFlag = true;

    elseif strcmpi(varargin{i},'IdsOnly')
        IdsOnly = true;
        SampleAnnotationFlag = false;
    elseif strcmpi(varargin{i},'Seperate')
        SeperateFiles = true;
    elseif strcmpi(varargin{i},'ProbeId')
        i = i + 1;
        ProbeIds = varargin{i};
        [~,~,ProbeIndx]=intersect(ProbeIds,DATA.ProbeId,'stable');
    elseif strcmpi(varargin{i},'SampleId')
        i = i + 1;
        SampleIds = varargin{i};
    end
end

format_str_txt = sprintf('%s%%s',Delimiter);
format_str_val = sprintf('%s%%g',Delimiter);
format_str_short = sprintf('%%s');

[fid,message] = fopen(FileOut,'w');
if  fid == -1
    disp(FileOut)
    disp(message)
    return
end

if SeperateFiles
    IdsOnly = true;
    [filepath,name,ext] = fileparts(FileOut);
    [fid_SA,message] = fopen(strcat(name,'_SampleAnnotation',ext),'w');
    if  fid_SA == -1
        disp(FileOut)
        disp(message)
        return
    end
    fprintf(fid_SA,'SampleIdId');
    if  ~isempty(DATA.RowAnnotationFields)
        fprintf(fid_SA,format_str_txt,DATA.RowAnnotationFields{:});
    end
    fprintf(fid_SA,'\n');
    for i=1:DATA.nRow
        fprintf(fid_SA,format_str_short,DATA.RowId{i});
        fprintf(fid_SA,format_str_txt,DATA.RowAnnotation{i,:});
        fprintf(fid_SA,'\n');
    end
    fclose(fid_SA);

    [fid_VA,message] = fopen(strcat(name,'_VariableAnnotation',ext),'w');
    if  fid_SA == -1
        disp(FileOut)
        disp(message)
        return
    end
    fprintf(fid_VA,'VariableId');
    if  ~isempty(DATA.ColAnnotationFields)
        fprintf(fid_VA,format_str_txt,DATA.ColAnnotationFields{:});
    end
    fprintf(fid_VA,'\n');
    for i=1:DATA.nCol
        fprintf(fid_VA,format_str_short,DATA.ColId{i});
        fprintf(fid_VA,format_str_txt,DATA.ColAnnotation{i,:});
        fprintf(fid_VA,'\n');
    end
    fclose(fid_VA);

end

fprintf(fid,'Id');
if SampleAnnotationFlag
    if  ~isempty(DATA.RowAnnotationFields)
        fprintf(fid,format_str_txt,DATA.RowAnnotationFields{:});
        fprintf(fid,'\t');
    end
end

if SurvivalFlag
    if  isfield(DATA,'SURVIVAL')
        for i=1:length(DATA.SURVIVAL.SurvivalTypes)
            fprintf(fid,'%s (Event)\t%s Time (%s)\t',DATA.SURVIVAL.SurvivalTypes{i},DATA.SURVIVAL.SurvivalTypes{i},DATA.SURVIVAL.Units{i});
        end
    end
end



if DATAFlag
    fprintf(fid,format_str_txt,DATA.ColId{:});
end

fprintf(fid,'\n');

for i=1:DATA.nRow
    fprintf(fid,format_str_short,DATA.RowId{i});
    if SampleAnnotationFlag
        if  ~isempty(DATA.RowAnnotation)
            fprintf(fid,format_str_txt,DATA.RowAnnotation{i,:});
            fprintf(fid,'\t');
        end
    end
    if SurvivalFlag
        for j=1:length(DATA.SURVIVAL.SurvivalTypes)
            fprintf(fid,'%s\t%g\t',DATA.SURVIVAL.SurvEvent{i,j},DATA.SURVIVAL.SurvTime(i,j));
        end
    end

    if DATAFlag

        fprintf(fid,format_str_val,DATA.X(i,:));
    end
    fprintf(fid,'\n');
end

fclose(fid);






end