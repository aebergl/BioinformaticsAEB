function fh = VolcanoPlotStacked(x_data,y_data,xlabel_txt,ylabel_txt)

FontSize=12;
min_alpha = 0.01;
max_alpha = 0.4;
min_size = 1;
max_size = 200;

fh=figure('Name','Volcano Plot','Color','w','Tag','Volcano Plot figure','Units','centimeters');
fh.Position(3:4) = [12 10];
ah = axes(fh,'NextPlot','add','tag','Volcano Plot','box','off','Layer','top','FontSize',FontSize,'YAxisLocation','origin');
ah.LineWidth = 74;
ah.XGrid= 'on';
ah.YGrid= 'on';

%uistack(ah,'top')


dist_val = pdist2([x_data y_data],[0 0],"seuclidean");
dist_alpha = rescale(dist_val.^2,min_alpha,max_alpha);
dist_size = rescale(dist_val.^2,min_size,max_size);

indx_pos = x_data > 0;
indx_neg = x_data < 0;

x_data_pos = x_data(indx_pos);
y_data_pos = y_data(indx_pos);
dist_pos_size = dist_size(indx_pos);
dist_pos_alpha = dist_alpha(indx_pos);

x_data_neg = x_data(indx_neg);
y_data_neg = y_data(indx_neg);
dist_neg_size = dist_size(indx_neg);
dist_neg_alpha = dist_alpha(indx_neg);

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

% min_x=min(x_data_neg,[],'all','omitnan')
% max_x=max(x_data_pos,[],'all','omitnan')
% nudge_x = (max_x-min_x)/20;

% ah.XLim=[min_x-nudge_x max_x+nudge_x];

