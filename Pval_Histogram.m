function fh = Pval_Histogram(P_values,nBins,Cmap,TitleTxt)
FontSize = 8;
Linewidth = 0.5;
FigureSize = [3,2];

fh = figure('Name','p-value histogram','Color','w','Tag','p-value histogram','Units','inches');
fh.Position(3:4) = FigureSize;
ah=axes(fh,'NextPlot','add','tag','p-value histogram','Box','off','FontSize',FontSize,'Linewidth',1);

hh=histogram(ah,P_values,nBins,'LineWidth',Linewidth,'Normalization','count','FaceColor',Cmap,...
    'FaceAlpha',0.9);
ah.XLabel.String = {'p-value'};
ah.XLabel.FontSize = FontSize;
ah.XLim = [-0.05 1.05];
ah.XTick = 0:0.1:1;
ah.XTickLabelRotation = -45;
ah.YLabel.String = 'Count';
ah.YLabel.FontSize = FontSize;
ah.YTick = [];

if ~isempty(TitleTxt)
   title(TitleTxt,'FontSize',FontSize,'FontWeight','Normal')
    
end

end
