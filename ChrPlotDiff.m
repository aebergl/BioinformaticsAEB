function fh = ChrPlotDiff(DATA,ChrPos,Chrfield,Chr,Y_Type,SizeType,ColorType,varargin)
% USAGE:
%   fh = ChrPlotDiff(DATA,ChrPos,Chrfield,Chr,Y_Type,SizeType,ColorType,varargin)
%   Add column annotation from file to DATA
%
% INPUTS:
% * DATA: DATA structure with results
% * ChrPos: Name of chr position field in RowAnnotationFields 'CpG_beg'
% * Chrfield: Name of Chr field in RowAnnotationFields 'CpG_chrm'
% * Chr: Name of Chromsome to be used 'chr1'
% * Y_Type: Name variable to be used on Y-axis from ColId 'HR coxreg DSS'
% * SizeType: Name variable to be used as size from ColId 'HR coxreg DSS'
% * ColorType: Name variable to be used as size from ColId 'p coxreg DSS'
%
% OUTPUTS:
% * fh: Figure handle to Chromosome figure
%   options ---------------------------------------
%
%   'GENES'         List of genes to mark and gene symbol field from RowAnnotationFields
%                   {'TMEM183A','RBBP5','ELK4'}, 'gene_HGNC'
%   'PlotRange'     Vector with start and stop position to to show in figure [203000000 226100000]
%   'MinMaxSize'    Vector with min and max marker size [1 80]
%   'FigSize'       Vector with figure width and hight in inches [6 2]
%   'REGION'        Highligth one or multiple regions cell structure with 
%                   {[Xstart Xstop], [Ystart Xstop], LabelTxt} 
%   'CytoBand'      Displays the Cytoband nased on HG38, give unit to be used, 'mb'
%   'YValCutOff'    Selects points with abs(value) larger than YValCutOff [0]
%   'SizeLegend'    Displays a size legend defined with the following input vectors:
%                   Size cut-off for the different sizes  [1 0.05 0.01 0.001]
%                   Marker size for the different cut-off [5 10 20 30]%
%                   Y position for the different markers sizes = [2.5 3 3.5 4]
%   'MarkSelected'  Highligh selected points
%   'Print'         Do not transpose the input file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Anders Berglund, 2020 aebergl at gmail.com                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Defaults
PlotRange = [];
FigSize = [6 2];
MinMaxSize = [1 80];
FontSize = 6;
RightMargin = 0.4;
TopMargin = 0;
printResults=false;
MarkSelected=false;
Delimiter = '\t';
pCeil = 15;
ChrUnit = 'mb';
YLabel = '';
Colorlabel = '';
SizeLable= '';
CytoBand = false;
GENES = [];
REGION = [];
YValCutOff = 0;
ColorValCutOff = 0;
CytoBandfile = 'cytoBandIdeo_38.txt';
SizeLegend = false;
% Check Input



i=0;
while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'GENES')
        i = i + 1;
        GeneField = varargin{i};
        i = i + 1;
        GENES = varargin{i};
    elseif strcmpi(varargin{i},'REGION')
        i = i + 1;
        REGION = varargin{i};           
    elseif strcmpi(varargin{i},'PlotRange')
        i = i + 1;
        PlotRange = varargin{i};          
    elseif strcmpi(varargin{i},'MinMaxSize')
        i = i + 1;
        MinMaxSize = varargin{i};
    elseif strcmpi(varargin{i},'FigSize')
        i = i + 1;
        FigSize = varargin{i};
    elseif strcmpi(varargin{i},'Print')
        printResults = true;
    elseif strcmpi(varargin{i},'CytoBand')
        CytoBand = true;
        i = i + 1;
        ChrUnit = varargin{i};
    elseif strcmpi(varargin{i},'YValCutOff')
        i = i + 1;
        YValCutOff = varargin{i};
    elseif strcmpi(varargin{i},'ColorValCutOff')
        i = i + 1;
        ColorValCutOff = varargin{i};
    elseif strcmpi(varargin{i},'MarkSelected')
        MarkSelected = true;
    elseif strcmpi(varargin{i},'SizeLegend')
        i = i + 1;
        SizeLegendCutOff = varargin{i};
        i = i + 1;
        LegendSizeVal = varargin{i};
        i = i + 1;
        LegendYPos = varargin{i};
        SizeLegend = true;
    end
end

ChrPosColumn = strcmpi(ChrPos,DATA.RowAnnotationFields);
if ~any(ChrPosColumn)
    error('Given value for ChrPos not found, %s not found in RowAnnotationFields',ChrPos{1})
end

ChrColumn = strcmpi(Chrfield,DATA.RowAnnotationFields);
if ~any(ChrColumn)
    error('Given value for Chrfield not found, %s not found in RowAnnotationFields',ChrColumn{1})
end
if ~isempty(GENES)
GeneColumn = strcmpi(GeneField,DATA.RowAnnotationFields);
if ~any(ChrColumn)
    error('Given value for GeneField not found, %s not found in RowAnnotationFields',ChrColumn{1})
end
end

%
% if DATA.NumProbes > 500000
%     ChipType = 'EPIC';
%     CytoBandfile = 'cytoBandIdeo_37.txt';
% else
%     ChipType='M450K';
%     CytoBandfile = 'cytoBandIdeo_37.txt';
% end


indx_Chr = strcmp(Chr,DATA.RowAnnotation(:,ChrColumn));
ChrPos = DATA.RowAnnotation(indx_Chr,ChrPosColumn);
FullAnnotation = DATA.RowAnnotation(indx_Chr,:);
Full_X_Data = DATA.X(indx_Chr,:);

% Remove if EPIC
try
    ChrPos = cellfun(@(x) str2num(x), ChrPos, 'UniformOutput', 0);
catch
end

ChrPos = cell2mat(ChrPos);
ChrPos = double(ChrPos);

switch ChrUnit
    case 'bp'
        ChrPos = ChrPos ./ 1.00;
        UnitVal = 1;
        UnitTxt = '(bp)';
    case 'kb'
        ChrPos = ChrPos ./ 1000.00;
        PlotRange = PlotRange ./ 1000.00;
        UnitVal = 2;
        UnitTxt = '(kb)';
    case 'mb'
        ChrPos = ChrPos ./ 1000000.00;
        PlotRange = PlotRange ./ 1000000.00;
        if ~isempty(REGION)
            REGION(:,1) = cellfun(@(x) x/1000000 ,REGION(:,1),'UniformOutput',false);
        end
        UnitVal = 3;
        UnitTxt = '(mb)';
end

% Select Y value
indx_Y_Val = strcmpi(Y_Type,DATA.ColId);
if any(indx_Y_Val)
    Y_Val = DATA.X(indx_Chr,indx_Y_Val);
else
    error('%s not found',Y_Type)
end

% Select Size value
indx_SizeVal = strcmpi(SizeType,DATA.ColId);
if any(indx_SizeVal)
    SizeVal = DATA.X(indx_Chr,indx_SizeVal);
else
    error('%s not found',SizeType)
end

% Select Color value
indx_ColorVal = strcmpi(ColorType,DATA.ColId);
if any(indx_SizeVal)
    ColorVal = DATA.X(indx_Chr,indx_ColorVal);
else
    error('%s not found',ColorType)
end

switch Y_Type
    case 'Delta Average'
        YLabel = {'\Delta \beta-value'};
    case 'HR coxreg OS'
        Y_Val =  log2(Y_Val);
        YLabel = {'log_2(HR OS)'};
    case 'HR coxreg DSS'
        Y_Val =  log2(Y_Val);
        YLabel = {'log_2(HR DSS)'};
    case 'HR coxreg PFI'
        Y_Val =  log2(Y_Val);
        YLabel = {'log_2(HR PFI)'};
    case 'r Spearman'
        YLabel = {'Spearman''s \rho'};
end

switch ColorType
    case 'q t-test'
        ColorVal =  -log10(ColorVal);
        Colorlabel = {'-log_1_0(q-value)'};
    case 'fdr t-test'
        ColorVal =  -log10(ColorVal);
        Colorlabel = {'-log_1_0(fdr p-value)'};
    case 'p coxreg DSS'
        ColorVal =  -log10(ColorVal);
        Colorlabel = {'-log_1_0(p coxreg DSS)'};
    case 'q coxreg OS'
        ColorVal =  -log10(ColorVal);
        Colorlabel = {'-log_1_0(q coxreg OS)'};

    case 'p coxreg PFI'
        ColorVal =  -log10(ColorVal);
        Colorlabel = {'-log_1_0(p coxreg PFI)'};
    case 'p Spearman'
        ColorVal =  -log10(ColorVal);
        Colorlabel = {'-log_1_0(p Spearman)'};
    case 'Range'
        Colorlabel = {'Range'};

end

switch SizeType
    case 'q t-test'
        SizeVal = -log10(SizeVal);
        SizeVal(SizeVal>pCeil) = pCeil;
    case 'Delta Average'
        SizeVal =  abs(SizeVal);
    case 'HR coxreg OS'
        SizeVal =  abs(log2(SizeVal));
    case 'HR coxreg DSS'
        SizeVal =  abs(log2(SizeVal));
    case 'HR coxreg PFI'
        SizeVal =  abs(log2(SizeVal));
    case 'p Spearman'
        SizeVal = -log10(SizeVal);

end
SizeVal =  abs(SizeVal);
% Select Range to plot. Do not work with CytoBand
if ~isempty(PlotRange)
    plot_indx = ChrPos >= PlotRange(1) & ChrPos<=PlotRange(2);
    ChrPos = ChrPos(plot_indx);
    Y_Val = Y_Val(plot_indx);
    SizeVal = SizeVal(plot_indx);
    ColorVal = ColorVal(plot_indx);
    FullAnnotation = FullAnnotation(plot_indx,:);
    Full_X_Data = Full_X_Data(plot_indx,:);
end
nGroups = length(Y_Val);
[~,sort_indx] = sort(ColorVal,'ascend');
ChrPos = ChrPos(sort_indx);
Y_Val = Y_Val(sort_indx);
SizeVal = SizeVal(sort_indx);
ColorVal = ColorVal(sort_indx);
FullAnnotation = FullAnnotation(sort_indx,:);
Full_X_Data = Full_X_Data(sort_indx,:);

SizeVal(SizeVal<0) = 0; % make sure there is no negative values
SizeValPlot = rescale(SizeVal,MinMaxSize(1),MinMaxSize(2));

% LegendSizeValMat = sum(repmat(SizeLegendCutOff,nGroups,1) >= repmat(SizeVal,1,length(SizeLegendCutOff)),2);
% SizeValPlot = LegendSizeVal(LegendSizeValMat);

indx_selected = abs(Y_Val) > YValCutOff & ColorVal > ColorValCutOff;

FullAnnotation_Selected = FullAnnotation(indx_selected,:);
Full_X_Data_Selected = Full_X_Data(indx_selected,:);
SizeValPlot_Selected = SizeValPlot(indx_selected,:);
Y_Val_Selected = Y_Val(indx_selected,:);
ChrPos_Selected = ChrPos(indx_selected,:);



fh = figure('Name','Chr. plot','Color','w','Tag','Chr. plot','GraphicsSmoothing','on','Unit','Inches');
fh.Position(3:4) = FigSize;
ah = axes(fh,'NextPlot','add','tag','Dot plot','Box','on','FontSize',FontSize,'Linewidth',0.5,'YGrid','on',...
    'Units','normalized','PositionConstraint','outerposition','Clipping','off');
cmap = colormap(colorcet('L17'));
colormap(cmap);

sh=scatter(ChrPos,Y_Val,SizeValPlot,ColorVal,'filled');
if MarkSelected
    sh_sel=scatter(ChrPos_Selected,Y_Val_Selected,SizeValPlot_Selected,[0 0 0]);
end

nudge_X = range(ah.XLim)/50;
nudge_Y = range(ah.YLim)/50;
ah.XLim = [min(ChrPos)-nudge_X max(ChrPos)+nudge_X];
ah.XAxis.TickDirection ='out';
ah.XAxis.TickLabelRotation = 0;
ah.XAxis.TickLength = [0.00500 0.0250];


for i=1:numel(GENES)
    gene=GENES{i};
    indx_1=strcmp(gene,FullAnnotation(:,GeneColumn));
    indx_2=contains(FullAnnotation(:,GeneColumn),[';',gene]);
    indx_3=contains(FullAnnotation(:,GeneColumn),[gene,';']);
    indx = indx_1 | indx_2 | indx_3;
    chr_pos_gene = ChrPos(indx);
    Y_Val_gene = Y_Val(indx);
    SizeValPlot_Gene = SizeValPlot(indx);
    y_line_val = max(Y_Val) + nudge_Y;
    y_line_val = max(Y_Val_gene) + nudge_Y;
    line([min(chr_pos_gene) max(chr_pos_gene)], [y_line_val y_line_val],'LineWidth',0.75,'Color','k','LineStyle','-');
    text((min(chr_pos_gene)+max(chr_pos_gene))/2,y_line_val ,gene,'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',FontSize-1,'FontAngle','italic')
    sh_sel=scatter(chr_pos_gene,Y_Val_gene,SizeValPlot_Gene,[0 0 0]);
end

for i=1:size(REGION,1)
    x = REGION{i,1};
    y = REGION{i,2};
    rectangle('Position',[x(1), y(1), x(2)-x(1), y(2)-y(1)],'LineWidth',0.5,'EdgeColor','k','LineStyle','-')
    text( (x(1) + x(2) )/2,y(2) ,REGION{i,3},'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',FontSize,'FontAngle','italic');
end

if CytoBand
    try
        chr_txt=strrep(Chr,'chr','');
        chromosomeplot(CytoBandfile, chr_txt, 'addtoplot', ah,'Orientation','Horizontal','unit', UnitVal);
        indx=findobj(fh.Children(2),'FontSize',8);
        [indx.FontSize]=deal(FontSize-2);
        ChrTxt = strcat('Chr. ', chr_txt);
        text(fh.Children(2),fh.Children(2).XLim(1)-0.005,0.7 ,ChrTxt,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',FontSize)
        TopMargin = 0.3;
    catch
    end
end
ylabel(ah,YLabel);

line(ah,ah.XLim,[0 0],'LineWidth',0.75,'Color','k','LineStyle','-');
ah.XAxis.TickLabels = strcat(ah.XAxis.TickLabels, UnitTxt);

ah.Units='inches';
ah.OuterPosition(3:4) = [FigSize(1)-RightMargin FigSize(2)-TopMargin];

ch = colorbar(ah,'Units','inches','FontSize',FontSize,...
    'Position',[ah.Position(1) + ah.Position(3)+0.1, ah.Position(2) 0.1, FigSize(2)/3]);
ch.Label.String=Colorlabel;
ch.FontSize=FontSize;



if SizeLegend

    shl = scatter(ah,ah.XLim(2)+nudge_X*1.5,LegendYPos,LegendSizeVal,[0 0 0],'filled');
    text(ah,ah.XLim(2)+nudge_X*1.5,ah.YLim(2),'p-value','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',FontSize)

    for i = 1:length(LegendSizeVal)
        if i==1
            txt_str = 'N.S.';
        else
            txt_str = num2str(SizeLegendCutOff(i));
        end
        text(ah,ah.XLim(2)+nudge_X*2.3,LegendYPos(i),txt_str,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',FontSize)
    end
end
if CytoBand
    fh.Children(2).Units='Inches';
    fh.Children(2).Position(3)=fh.Children(4).Position(3);
end
if printResults
    format_str_txt = sprintf('%s%%s',Delimiter);
    format_str_val = sprintf('%s%%g',Delimiter);
    %format_str_short = sprintf('%%s');
    fprintf(format_str_txt,DATA.RowAnnotationFields{:});
    fprintf(format_str_txt,DATA.ColId{:});
    fprintf('\n');
    for i=1:size(FullAnnotation_Selected,1)
        fprintf(format_str_txt,FullAnnotation_Selected{i,:});
        fprintf(format_str_val,Full_X_Data_Selected(i,:));
        fprintf('\n');
    end
end


