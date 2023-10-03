function fh = PlotDensityHistogramDATA(DATA,Type,CMap)

LineWidth = 1;
FontSize = 12;
nBins = 100;
BandwidthValue = 0.05;
nPoints=10;
KernalDensity = false;
AlphaValue  = 0.2;

switch lower(Type)
    case 'samples'
        X = DATA.X;

    case 'variables'
        X = DATA.X';
end
[nRow,nCol] = size(X);

Xmean = mean(X,2,"omitnan");
[~,indx] = sort(Xmean,'ascend');

X = X(indx,:);

N=zeros(nRow,nBins);


for i=1:nRow
    %[N(i,:),edges] = ksdensity(X(i,:),'npoints',nBins,'Bandwidth',BandwidthValue);
    [N(i,:),edges] = histcounts(X(i,:),nBins);

end

fh = figure('Name','Density Histogram Plot','Color','w','Tag','Density Histogram Plot','Units','inches','Colormap',CMap);

ah = axes(fh,'NextPlot','add','tag','Density Histogram Plot','Box','on','FontSize',FontSize,'Linewidth',0.5,...
    'ActivePositionProperty','outerposition','XGrid','on','YGrid','on');

N=normalize(N,2,"range");

[px,py] = GetPatch(nRow,nBins);

patch(ah, px/nBins, py, reshape(N,1,nRow*nBins), 'linestyle', 'none');
colormap(ah,CMap);

% caxis(ax_HM,options.CLim);
ah.XLim = [.5/nBins (nBins+.5)/nBins];
ah.YLim = [.5 nRow+.5];
% 
% fcb_HM=colorbar('peer',ax_HM,'location','EastOutside','Tag','Colorbar for HeatMap');
% fcb_HM.Position=[0.92 0.5 0.03 0.25];
% fcb_HM.Limits = options.CLim;
% fcb_HM.Ticks = options.CLim(1):options.CLim(2);
% fcb_HM.FontSize = 14;