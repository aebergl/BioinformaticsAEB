function DATA = Geo2DATA(GeoId,SampleFile,ProbeAnnotationFile)

try
    fprintf('Loading GEO dataset\n')
    GEO_DATA = getgeodata(GeoId);
catch
    error('Could not download GEO datset: %s\nPlease check name or connection\n',GeoId)
end

DATA.X = (GEO_DATA.Data.(':')(':'))';
[nSamples, nVar] = size(DATA.X);

DATA.NumSamples = nSamples;
DATA.NumProbes = nVar;
DATA.ProbeId = strtrim(GEO_DATA.Data.RowNames);
DATA.SampleId = strtrim(GEO_DATA.Data.ColNames);

% Attemting to load Probe Annotation file
PlatformId = GEO_DATA.Header.Series.platform_id;
if ~isempty(PlatformId)
    fprintf('Loading platform annotation: %s\n',PlatformId)
    try
        GPL_DATA = getgeodata(PlatformId);
    catch
        error('Could not download PlatformId data: %s\nPlease check name and connection\n',PlatformId)
    end
    ProbeAnnotation = cell(DATA.NumProbes,numel(GPL_DATA.ColumnNames));
    ProbeAnnotation(:) = {'---'};
    
    AnnotationTxt = GPL_DATA.Data;
    AnnotationTxt = cellfun(@(x) num2str(x),AnnotationTxt,'UniformOutput',false);
    [~,indx1,indx2] = intersect(DATA.ProbeId,AnnotationTxt(:,1),'Stable');
    ProbeAnnotation(indx1,:) = AnnotationTxt(indx2,:);
    if strcmpi('GPL5175',PlatformId)
        fprintf('Parsing platform annotation: %s\n',PlatformId)
        [DATA.ProbeAnnotation, DATA.ProbeAnnotationFields] = ParseGPL5175(ProbeAnnotation,GPL_DATA.ColumnNames);
    else
        DATA.ProbeAnnotationFields = GPL_DATA.ColumnNames;
        DATA.ProbeAnnotation = ProbeAnnotation;
    end
end
[SampleInfo, ColumnNames] = Geo2SampleInfo(GEO_DATA,SampleFile);

DATA.SampleAnnotationFields = ColumnNames;
DATA.SampleAnnotation = SampleInfo;

DATA.ProbeXDim = 2;
DATA.SampleXDim = 1;
DATA.Type = 'Data';
