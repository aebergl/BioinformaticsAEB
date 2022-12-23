function fh = ChrPlotDiff(DATA,Chr,PlotRange,Y_Type,SizeType,ColorType,varargin)
FontSize = 10;
minSize = 1;
maxSize = 150;
printResults=false;
MarkSelected=false;
Delimiter = '\t';
pCeil = 15;
ChrUnit = 'bp';
% ChrUnit = 'kb';
ChrUnit = 'mb';
YLabel = '';
Colorlabel = '';
SizeLable= '';
CytoBand = false;
GENES=[];
REGION = [];
YValCutOff = 0;
ColorValCutOff = 0;
% Check Input
i=0;
while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'GENES')
        i = i + 1;
        GENES = varargin{i};
    elseif strcmpi(varargin{i},'REGION')
        i = i + 1;
        REGION = varargin{i};
    elseif strcmpi(varargin{i},'Print')
        printResults = true;
    elseif strcmpi(varargin{i},'ChrUnit')
        i = i + 1;
        ChrUnit = varargin{i};
    elseif strcmpi(varargin{i},'CytoBand')
        CytoBand = true;
    elseif strcmpi(varargin{i},'YValCutOff')
        i = i + 1;
        YValCutOff = varargin{i};
    elseif strcmpi(varargin{i},'ColorValCutOff')
        i = i + 1;
        ColorValCutOff = varargin{i};
    elseif strcmpi(varargin{i},'MarkSelected')
        MarkSelected = true;
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
CytoBandfile = 'cytoBandIdeo_38.txt';
ChrColumn = strcmpi('CHR',DATA.RowAnnotationFields);
ChrPosColumn = strcmpi('MAPINFO',DATA.RowAnnotationFields);
ChrColumn = strcmpi('CpG_chrm',DATA.RowAnnotationFields);
ChrPosColumn = strcmpi('CpG_beg',DATA.RowAnnotationFields);
GeneColumn = strcmpi('gene_HGNC',DATA.RowAnnotationFields);
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
    case 'HR coxreg DSS'
        Y_Val =  log2(Y_Val);
        YLabel = {'log_2(HR DSS)'};
    case 'HR coxreg PFI'
        Y_Val =  log2(Y_Val);
        YLabel = {'log_2(HR PFI)'};
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
    case 'p coxreg PFI'
        ColorVal =  -log10(ColorVal);
        Colorlabel = {'-log_1_0(p coxreg PFI)'};        
end

switch SizeType
    case 'q t-test'
        SizeVal = -log10(DATA.RES(nRES).pTT(indx_Chr));
        SizeVal(SizeVal>pCeil) = pCeil;
    case 'Delta Average'
        SizeVal =  abs(SizeVal);
    case 'HR coxreg DSS'
        SizeVal =  abs(log2(SizeVal));
      case 'HR coxreg PFI'
        SizeVal =  abs(log2(SizeVal));
end
SizeVal =  abs(SizeVal);
if ~isempty(PlotRange)
    plot_indx = ChrPos >= PlotRange(1) & ChrPos<=PlotRange(2);
    ChrPos = ChrPos(plot_indx);
    Y_Val = Y_Val(plot_indx);
    SizeVal = SizeVal(plot_indx);
    ColorVal = ColorVal(plot_indx);
    FullAnnotation = FullAnnotation(plot_indx,:);
    Full_X_Data = Full_X_Data(plot_indx,:);
end

[~,sort_indx] = sort(ColorVal,'ascend');
ChrPos = ChrPos(sort_indx);
Y_Val = Y_Val(sort_indx);
SizeVal = SizeVal(sort_indx);
ColorVal = ColorVal(sort_indx);
FullAnnotation = FullAnnotation(sort_indx,:);
Full_X_Data = Full_X_Data(sort_indx,:);

SizeVal(SizeVal<0) = 0; % make sure there is no negative values
SizeValPlot = rescale(SizeVal,minSize,maxSize);

indx_selected = abs(Y_Val) > YValCutOff & ColorVal > ColorValCutOff;

FullAnnotation_Selected = FullAnnotation(indx_selected,:);
Full_X_Data_Selected = Full_X_Data(indx_selected,:);
SizeValPlot_Selected = SizeValPlot(indx_selected,:);
Y_Val_Selected = Y_Val(indx_selected,:);
ChrPos_Selected = ChrPos(indx_selected,:);



fh = figure('Name','Chr. plot','Color','w','Tag','Chr. plot','GraphicsSmoothing','on');
fh.Position(3:4) = [700 250];
ah = axes(fh,'NextPlot','add','tag','Dot plot','Box','on','FontSize',FontSize,'Linewidth',0.75,'YGrid','on');
%ah.Position = [0.07 0.17 0.8 0.9];
cmap = colormap(colorcet('L17'));

%cmap = flipud(cmap);
colormap(cmap);
sh=scatter(ChrPos,Y_Val,SizeValPlot,ColorVal,'filled');
if MarkSelected
sh_sel=scatter(ChrPos_Selected,Y_Val_Selected,SizeValPlot_Selected,[0 0 0]);
end
ch = colorbar(ah);
ch.Label.String=Colorlabel;
nudge_X = range(ah.XLim)/100;
nudge_Y = range(ah.YLim)/100;
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
    y_line_val = max(Y_Val) + nudge_Y;
    line([min(chr_pos_gene) max(chr_pos_gene)], [y_line_val y_line_val],'LineWidth',0.75,'Color','k','LineStyle','-');
    text((min(chr_pos_gene)+max(chr_pos_gene))/2,y_line_val ,gene,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',FontSize,'FontAngle','italic')
end

for i=1:size(REGION,1)
    line([REGION{i,1}(1) REGION{i,1}(2)], [REGION{i,2}(1) REGION{i,2}(1)],'LineWidth',0.75,'Color','k','LineStyle','-');
    text( (REGION{i,1}(1) + REGION{i,1}(2) )/2,REGION{i,2} ,REGION{i,3},'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',FontSize,'FontAngle','italic');
end

if CytoBand
    try
        chr_txt=strrep(Chr,'chr','');
        chromosomeplot(CytoBandfile, chr_txt, 'addtoplot', ah,'Orientation','Horizontal','unit', UnitVal);
        indx=findobj(fh.Children(2),'FontSize',8);
        [indx.FontSize]=deal(7);
        ChrTxt = strcat('Chr. ', chr_txt);
        text(fh.Children(2),fh.Children(2).XLim(1)-0.005,0.7 ,ChrTxt,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',FontSize)

    catch
    end
end
ylabel(ah,YLabel);
%title(ah,sprintf('Chromosome %s',Chr),'FontSize',12,'FontWeight','Normal');
line(ah,ah.XLim,[0 0],'LineWidth',0.75,'Color','k','LineStyle','-');
ah.XAxis.TickLabels = strcat(ah.XAxis.TickLabels, UnitTxt);



if printResults
    format_str_txt = sprintf('%s%%s',Delimiter);
    format_str_val = sprintf('%s%%g',Delimiter);
    format_str_short = sprintf('%%s');
    fprintf(format_str_txt,DATA.RowAnnotationFields{:});
    fprintf(format_str_txt,DATA.ColId{:});
    fprintf('\n');
    for i=1:size(FullAnnotation_Selected,1)
        fprintf(format_str_txt,FullAnnotation_Selected{i,:});
        fprintf(format_str_val,Full_X_Data_Selected(i,:));
        fprintf('\n');
    end
end



