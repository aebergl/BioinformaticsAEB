function fh = VolcanoPlotStacked(x_data,y_data,xlabel_txt,ylabel_txt)

FontSize=12;
max_alpha = 0.4;
max_size = 200;

fh=figure('Name','Volcano Plot','Color','w','Tag','Volcano Plot figure','Units','centimeters');
fh.Position(3:4) = [12 10];
ah = axes(fh,'NextPlot','add','tag','Volcano Plot','box','on','Layer','top','FontSize',FontSize);
ah.LineWidth = 1;
ah.XGrid= 'on';
ah.YGrid= 'on';

uistack(ah,'top')

indx_pos = x_data > 0;
indx_neg = x_data < 0;

x_data_pos = x_data(indx_pos);
y_data_pos = y_data(indx_pos);

x_data_neg = -x_data(indx_neg);
y_data_neg = y_data(indx_neg);

dist_pos = pdist2([x_data_pos y_data_pos],[0 0],"seuclidean");
dist_pos_alpha = rescale(dist_pos.^2,0.01,max_alpha);
dist_pos_size = rescale(dist_pos.^2,1,max_size);

dist_neg  = pdist2([x_data_neg y_data_neg],[0 0],"seuclidean");
dist_neg_alpha = rescale(dist_neg.^2,0.01,max_alpha);
dist_neg_size = rescale(dist_neg.^2,1,max_size);

sh_pos = scatter(x_data_pos,y_data_pos,dist_pos_size,[1 0 0],'filled','Parent',ah);
sh_pos.AlphaDataMapping = 'none';
sh_pos.AlphaData = dist_pos_alpha;
sh_pos.MarkerFaceAlpha = 'flat';

sh_neg = scatter(x_data_neg,y_data_neg,dist_neg_size,[0 0 1],'filled','Parent',ah);
sh_neg.AlphaDataMapping = 'none';
sh_neg.AlphaData = dist_neg_alpha;
sh_neg.MarkerFaceAlpha = 'flat';


xlabel(ah,xlabel_txt);
ylabel(ah,ylabel_txt);
