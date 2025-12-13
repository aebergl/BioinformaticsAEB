function fh = Pval_Histogram(P_values,nBins,Cmap,TitleTxt)
fh = figure('Name','p-value histogram','Color','w','Tag','p-value histogram','GraphicsSmoothing','off');
ah=axes(fh,'NextPlot','add','tag','p-value histogram','Box','off','FontSize',16,'Linewidth',1);

hh=histogram(ah,P_values,nBins,'LineWidth',1,'Normalization','count','FaceColor',Cmap);
ah.XLabel.String = {'p-value'};
ah.XLabel.FontSize = 18;
ah.XLim = [-0.05 1.05];
ah.XTick = 0:0.1:1;

ah.YLabel.String = 'Count';
ah.YLabel.FontSize = 18;
ah.YTick = [];

if ~isempty(TitleTxt)
   title(TitleTxt,'FontSize',20,'FontWeight','Normal')
    
end

end
