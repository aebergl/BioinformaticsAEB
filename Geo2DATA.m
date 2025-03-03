function DATA = Geo2DATA(GeoId,SampleFile)
% USAGE:
%   fh  = PlotBoxPlotDATA(DATA,VariableId,GroupVariableName,GroupsToUse,varargin)
%   creates a scatter boxplot plot from DATA structure
%
% INPUTS:
% * DATA:       DATA structure w
% * VariableId: Id to be used for groups along the X-axis from RowAnnotationFields.
% * GroupVariableName:    List of Groups to be used and order for the X-axis, [] uses all
% * GroupsToUse:    Id to be used for groups along the Y-axis from RowAnnotationFields.
%
% OUTPUTS:
%   fh: Figure handle to dot plot figure
%
%   options ---------------------------------------
%
%   'FontSize'          FontSize for all text in figure [7]
%   'FigSize'           Vector with figure width and hight in inches [ [2 3.5] ]
%   'MarkerSize'        Marker size for scatter points [30]
%   'MarkerLineWidth'   Marker line width for scatter points [1]
%   'MarkerType'        Marker types for scatter points,{'o','d','^'} [{'o'}]
%   'MarkerFaceColor'   Marker edge color for scatter points ['none']
%   'ColorMap'          Marker edge color for scatter points [GetPalette('Science')]
%   'XJitterWidth'      Width for scatter points [0.6]
%   'BoxLineWidth'      LineWidth for boxes [1]
%   'BoxColor'          Color for box lines [ [0 0 0] ]
%   'BoxWidth'          Width for boxes [0.8]
%   'YlabelTxt'         Text to be used for Y label, if not defined VariableId will be used
%   'TitleTxt'          text for Title
%   'CalcStats'         Calculate stats between groups, Nx2 matrix defines comparisons, [] all
%   'StatType'          Type of comparison, Mann Whitney (MW) or T-test (t-test) ['MW']
%   'PlotStars'         Use starts instead of p-values
%   'StatLineWidth'     Line witth for line btween groups [0.5]
%   'Show_NS'           Show not significant results [false]
%   'TargetAxes'        Axes handle for target axist []
%   'XTickAngle'        Angle for X ticks [-45]
%   'SortData'          Sort group order based on median value
%   'MultipleY'         How to summarize multiple y's Mean' or 'PCA' ['Mean']
%   'DataTipId'         Id to be used for datatip 'RowId' or id from 'RowAnnotationFields' ['RowId']

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Anders Berglund, 2025 aebergl at gmail.com                           %

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
