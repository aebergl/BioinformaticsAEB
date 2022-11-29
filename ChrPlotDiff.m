function fh = ChrPlotDiff(DATA,Chr,PlotRange,Y_Type,SizeType,ColorType)
FontSize = 10;
minSize = 10;
maxSize = 150;
printResults=false;
Delimiter = '\t';

GENES={'MEG3','MEG8','DLK1','MIR376B'};

%ChrUnit = 'bp';
% ChrUnit = 'kb';
ChrUnit = 'mb';
%PlotRange = [];
%PlotRange = [57000000 59000000];
pCeil = 15;


YLabel = '';
Colorlabel = '';
SizeLable= '';
% 
% if DATA.NumProbes > 500000
%     ChipType = 'EPIC';
%     CytoBandfile = 'cytoBandIdeo_37.txt';
% else
%     ChipType='M450K';
%     CytoBandfile = 'cytoBandIdeo_37.txt';
% end
CytoBandfile = 'cytoBandIdeo_38.txt';
% ChrColumn = strcmpi('CHR',DATA.ProbeAnnotationColumns);
% ChrPosColumn = strcmpi('MAPINFO',DATA.ProbeAnnotationColumns);
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
        UnitVal = 2;
        UnitTxt = '(kb)';
    case 'mb'
        ChrPos = ChrPos ./ 1000000.00;
        PlotRange = PlotRange ./ 1000000.00;
        UnitVal = 3;
        UnitTxt = ' (mb)';
end

% Select Y value
indx_Y_Val = strcmpi(Y_Type,DATA.ColId);
if any(indx_Y_Val)
    Y_Val = DATA.X(indx_Chr,indx_Y_Val);
else
    
end

% Select Size value
indx_SizeVal = strcmpi(SizeType,DATA.ColId);
if any(indx_SizeVal)
    SizeVal = DATA.X(indx_Chr,indx_SizeVal);
else
    
end

% Select Color value
indx_ColorVal = strcmpi(ColorType,DATA.ColId);
if any(indx_SizeVal)
    ColorVal = DATA.X(indx_Chr,indx_ColorVal);
else
    
end

switch Y_Type
    case 'p t-test'
        Y_Val = -log10(DATA.RES(nRES).pTT(indx_Chr));
    case 'q t-test'
        Y_Val = -log10(DATA.RES(nRES).Q_TT(indx_Chr));
    case 'r_Spearman'
        Y_Val = (DATA.RES(nRES).r_Spearman(indx_Chr));
    case 'r_Pearson'
        Y_Val = (DATA.RES(nRES).r_Pearson(indx_Chr));
    case 'Range'
        Y_Val = DATA.RES(nRES).Range(indx_Chr);
    case 'p_logrank'
        Y_Val = DATA.RES(nRES).p_logrank(indx_Chr);
    case 'Delta Average'
        YLabel = {'\Delta \beta-value'};
end

switch ColorType
    case 'q t-test'
       ColorVal =  -log10(ColorVal);
       Colorlabel = {'-log_1_0(q-value)'};
    case 'fdr t-test'
       ColorVal =  -log10(ColorVal);
       Colorlabel = {'-log_1_0(fdr p-value)'};


    case 'p t-test'
        ColorVal = -log10(DATA.RES(nRES).Q_TT(indx_Chr));
        Colorlabel = {'-log_1_0(p-value)'};
    case 'r_Spearman'
        ColorVal = (DATA.RES(nRES).r_Spearman(indx_Chr)).^2;
    case 'r_Pearson'
        ColorVal = (DATA.RES(nRES).r_Pearson(indx_Chr)).^2;
    case 'Range'
        ColorVal = DATA.RES(nRES).Range(indx_Chr);
     case 'p_logrank'
        ColorVal = DATA.RES(nRES).p_logrank(indx_Chr);
     case 'MeanDiff'
        ColorVal = (DATA.RES(nRES).MeanDiff(indx_Chr));
     case 'p_Spearman'
        ColorVal = (DATA.RES(nRES).p_Spearman(indx_Chr));
end

switch SizeType
    case 'q t-test'
        SizeVal = -log10(DATA.RES(nRES).pTT(indx_Chr));
        SizeVal(SizeVal>pCeil) = pCeil;
    case 'p t-test'
        SizeVal = (DATA.RES(nRES).Q_TT(indx_Chr));
    case 'r_Spearman'
        SizeVal = (DATA.RES(nRES).r_Spearman(indx_Chr)).^2;
    case 'r_Pearson'
        SizeVal = (DATA.RES(nRES).r_Pearson(indx_Chr)).^2;
    case 'Range'
        SizeVal = DATA.RES(nRES).Range(indx_Chr);
     case 'p_logrank'
        SizeVal = DATA.RES(nRES).p_logrank(indx_Chr);
    case 'Delta Average'
        SizeVal =  abs(SizeVal);
     case 'p_Spearman'
        SizeVal = (DATA.RES(nRES).p_Spearman(indx_Chr));
end


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

fh = figure('Name','Chr. plot','Color','w','Tag','Chr. plot','GraphicsSmoothing','on');
fh.Position(3:4) = [1000 430];
ah = axes(fh,'NextPlot','add','tag','Dot plot','Box','on','FontSize',FontSize,'Linewidth',0.75,'YGrid','on');
%ah.Position = [0.07 0.17 0.8 0.9];
cmap = colormap(colorcet('L17'));

%cmap = flipud(cmap);
colormap(cmap);
sh=scatter(ChrPos,Y_Val,SizeValPlot,ColorVal,'filled');
ch = colorbar(ah);
ch.Label.String=Colorlabel;
nudge_X = range(ah.XLim)/100;
nudge_Y = range(ah.YLim)/100;
ah.XLim = [min(ChrPos)-nudge_X max(ChrPos)+nudge_X];
ah.XAxis.TickDirection ='out';
ah.XAxis.TickLabelRotation = 0;
ah.XAxis.TickLength = [0.00500 0.0250];

% try
% chromosomeplot(CytoBandfile, '14', 'addtoplot', ah,'Orientation','Horizontal','unit', UnitVal);
% indx=findobj(fh.Children(2),'FontSize',8);
% [indx.FontSize]=deal(7);
% ChrTxt = strcat('Chr. ', Chr);
% text(fh.Children(2),fh.Children(2).XLim(1)-0.005,0.7 ,ChrTxt,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',FontSize)
% 
% catch
% end

ylabel(ah,YLabel);
%title(ah,sprintf('Chromosome %s',Chr),'FontSize',12,'FontWeight','Normal');
line(ah,ah.XLim,[0 0],'LineWidth',0.75,'Color','k','LineStyle','-');
ah.XAxis.TickLabels = strcat(ah.XAxis.TickLabels, UnitTxt);

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


if printResults
    format_str_txt = sprintf('%s%%s',Delimiter);
    format_str_val = sprintf('%s%%g',Delimiter);
    format_str_short = sprintf('%%s');
    fprintf(format_str_txt,DATA.RowAnnotationFields{:});
    fprintf(format_str_txt,DATA.ColId{:});
    fprintf('\n');
    for i=1:size(FullAnnotation,1)
        fprintf(format_str_txt,FullAnnotation{i,:});
        fprintf(format_str_val,Full_X_Data(i,:));
        fprintf('\n');
    end


end



