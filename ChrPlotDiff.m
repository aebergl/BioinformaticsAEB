function fh = ChrPlotDiff(DATA,Chr,PlotRange,Y_Type,SizeType,ColorType)
FontSize = 10;
FontSize = 12;
minSize = 10;
maxSize = 150;
ChrUnit = 'bp';
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
indx = strcmp(Chr,DATA.RowAnnotation(:,ChrColumn));
ChrPos = DATA.RowAnnotation(indx,ChrPosColumn);

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
    Y_Val = DATA.X(indx,indx_Y_Val);
else
    
end

% Select Size value
indx_SizeVal = strcmpi(SizeType,DATA.ColId);
if any(indx_SizeVal)
    SizeVal = DATA.X(indx,indx_SizeVal);
else
    
end

% Select Color value
indx_ColorVal = strcmpi(ColorType,DATA.ColId);
if any(indx_SizeVal)
    ColorVal = DATA.X(indx,indx_ColorVal);
else
    
end

switch Y_Type
    case 'p t-test'
        Y_Val = -log10(DATA.RES(nRES).pTT(indx));
        Y_Val(Y_Val>pCeil) = pCeil;
    case 'q t-test'
        Y_Val = -log10(DATA.RES(nRES).Q_TT(indx));
    case 'r_Spearman'
        Y_Val = (DATA.RES(nRES).r_Spearman(indx));
    case 'r_Pearson'
        Y_Val = (DATA.RES(nRES).r_Pearson(indx));
    case 'Range'
        Y_Val = DATA.RES(nRES).Range(indx);
    case 'p_logrank'
        Y_Val = DATA.RES(nRES).p_logrank(indx);
    case 'Delta Average'
        YLabel = {'\Delta \beta-value'};
end

switch ColorType
    case 'q t-test'
        ColorVal =  -log10(ColorVal);
       Colorlabel = {'-log_1_0(q-value)'};
    case 'p t-test'
        ColorVal = -log10(DATA.RES(nRES).Q_TT(indx));
        Colorlabel = {'-log_1_0(p-value)'};
    case 'r_Spearman'
        ColorVal = (DATA.RES(nRES).r_Spearman(indx)).^2;
    case 'r_Pearson'
        ColorVal = (DATA.RES(nRES).r_Pearson(indx)).^2;
    case 'Range'
        ColorVal = DATA.RES(nRES).Range(indx);
     case 'p_logrank'
        ColorVal = DATA.RES(nRES).p_logrank(indx);
     case 'MeanDiff'
        ColorVal = (DATA.RES(nRES).MeanDiff(indx));
     case 'p_Spearman'
        ColorVal = (DATA.RES(nRES).p_Spearman(indx));
end

switch SizeType
    case 'q t-test'
        SizeVal = -log10(DATA.RES(nRES).pTT(indx));
        SizeVal(SizeVal>pCeil) = pCeil;
    case 'p t-test'
        SizeVal = (DATA.RES(nRES).Q_TT(indx));
    case 'r_Spearman'
        SizeVal = (DATA.RES(nRES).r_Spearman(indx)).^2;
    case 'r_Pearson'
        SizeVal = (DATA.RES(nRES).r_Pearson(indx)).^2;
    case 'Range'
        SizeVal = DATA.RES(nRES).Range(indx);
     case 'p_logrank'
        SizeVal = DATA.RES(nRES).p_logrank(indx);
    case 'Delta Average'
        SizeVal =  abs(SizeVal);
     case 'p_Spearman'
        SizeVal = (DATA.RES(nRES).p_Spearman(indx));
end

[~,sort_indx] = sort(ColorVal,'ascend');
ChrPos = ChrPos(sort_indx);
Y_Val = Y_Val(sort_indx);
SizeVal = SizeVal(sort_indx);
ColorVal = ColorVal(sort_indx);

if ~isempty(PlotRange)
    plot_indx = ChrPos >= PlotRange(1) & ChrPos<=PlotRange(2);
    ChrPos = ChrPos(plot_indx);
    Y_Val = Y_Val(plot_indx);
    SizeVal = SizeVal(plot_indx);
    ColorVal = ColorVal(plot_indx);
end

SizeVal(SizeVal<0) = 0; % make sure there is no negative values
SizeVal = rescale(SizeVal,minSize,maxSize);

fh = figure('Name','Chr. plot','Color','w','Tag','Chr. plot','GraphicsSmoothing','on');
fh.Position(3:4) = [1000 430];
ah = axes(fh,'NextPlot','add','tag','Dot plot','Box','on','FontSize',FontSize,'Linewidth',0.75,'YGrid','on');
%ah.Position = [0.07 0.17 0.8 0.9];
cmap = colormap(colorcet('L17'));

%cmap = flipud(cmap);
colormap(cmap);
sh=scatter(ChrPos,Y_Val,SizeVal,ColorVal,'filled');
ch = colorbar(ah);
ch.Label.String=Colorlabel;
nudge_X = range(ah.XLim)/100;
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


