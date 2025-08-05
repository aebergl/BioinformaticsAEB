function [] = PanCan_Survival_All(DATA,VariableId,varargin)
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

VariableIdentifier = false;
TumorTypesToUse = [];
TumorTypeColumnId = 'TCGA_Code';
SurvivalTypes = [];
TissueColumnId = "Sample_Type";
CutPoint = 'median';

% Check input
if nargin > 2
    ArgsList = {'VariableIdentifier','TumorTypesToUse','SurvivalTypes','CutPoint'};
    for j=1:2:numel(varargin)

        ArgType = varargin{j};
        ArgVal = varargin{j+1};
        if ~strncmpi(ArgType,ArgsList,numel(ArgType))
            error('Invalid input option: %s', ArgType);
        else
            switch lower(ArgType)
                case 'variableidentifier'
                    VariableIdentifier = ArgVal;
                case 'tumortypestouse'
                    TumorTypesToUse = ArgVal;
                case 'survivaltypes'
                    SurvivalTypes = ArgVal;
               case 'cutpoint'
                    CutPoint = ArgVal;

            end
        end
    end
end

% Get Tumor Types To Use
indx_TumorType = strcmpi(TumorTypeColumnId,DATA.RowAnnotationFields);
if ~any(indx_TumorType)
    error('Error. \n%s not found in DATA.RowAnnotationFields',TumorTypeColumnId);
elseif sum(indx_TumorType) > 1
    error('Warning. \nMultiple matches for %s found in DATA.RowAnnotationFields',TumorTypeColumnId);
else
    if isempty(TumorTypesToUse)
        TumorTypesToUse = unique(DATA.RowAnnotation(:,indx_TumorType));
    end
end

DATA = EditSamplesDATA(DATA,TumorTypesToUse,'Keep','SampleIdentifier',TumorTypeColumnId);
nTumorTypes = length(TumorTypesToUse);

% Selection of Y variable
if VariableIdentifier
    indx_VarId = strcmpi(VariableIdentifier,DATA.ColAnnotationFields);
    if ~any(indx_VarId)
        error('Error. \n%s not found in DATA.ColAnnotationFields',VariableIdentifier);
    elseif sum(indx_VarId) > 1
        error('Warning. \nMultiple matches for %s found in DATA.ColAnnotationFields',VariableIdentifier);
    else
        DATA = EditVariablesDATA(DATA,VariableId,'Keep','VariableIdentifier',VariableIdentifier);
    end
else
    DATA = EditVariablesDATA(DATA,VariableId,'Keep');
end

if isempty(SurvivalTypes)
    SurvivalTypes = DATA.SURVIVAL.SurvivalTypes;
end
nSurvivalTypes = length(SurvivalTypes);

CutPoint
fprintf('\n')
fprintf('Tumor Type')
fprintf('\t%s',SurvivalTypes{:})
fprintf('\n')
for i = 1:nTumorTypes
    fprintf('%s',TumorTypesToUse{i})
    DATA_Calc = EditSamplesDATA(DATA,TumorTypesToUse(i),'Keep','SampleIdentifier',TumorTypeColumnId);

    switch lower(TumorTypesToUse{i})
        case 'skcm'
            DATA_Calc = EditSamplesDATA(DATA_Calc,{'Primary solid Tumor','Metastatic'},'Keep','SampleIdentifier',TissueColumnId);
        case 'laml'
            DATA_Calc = EditSamplesDATA(DATA_Calc,{'Primary Blood Derived Cancer - Peripheral Blood'},'Keep','SampleIdentifier',TissueColumnId);
        otherwise
            DATA_Calc = EditSamplesDATA(DATA_Calc,{'Primary solid Tumor'},'Keep','SampleIdentifier',TissueColumnId);
    end

    for j = 1:nSurvivalTypes

        indx_SurvivalType = strcmp(SurvivalTypes(j),DATA_Calc.SURVIVAL.SurvivalTypes);
        TimeVar = DATA_Calc.SURVIVAL.SurvTime(:,indx_SurvivalType);
        EventVar = DATA_Calc.SURVIVAL.SurvEvent(:,indx_SurvivalType);
        switch lower(DATA_Calc.Info.Type)
            case 'gistic'
                x = DATA_Calc.X;
                x = num2cell(x);
                x = cellfun(@(x) num2str(x), x,'UniformOutput',false);
                try 
                [p,fH,stats] = MatSurv( TimeVar, EventVar,x,'Timeunit','Months','NoPlot',1,'NoWarnings',1,'Print',0);
                catch 
                    p=NaN;
                end
                fprintf('\t%g',p)
            case 'geneexpression'
                 x = DATA_Calc.X;
                 try 
                [p,fH,stats] = MatSurv( TimeVar, EventVar,x,'Timeunit','Months','NoPlot',1,'NoWarnings',1,'Print',0,'CutPoint',CutPoint);
                 catch
                     p=NaN;
                 end
                fprintf('\t%g',p)
            otherwise
        end
        
    end
fprintf('\n')
end
