function fh = VolcanoPlotStacked(x_data,y_data,xlabel_txt,ylabel_txt)

PlotType = 'simple';
x_CutOff = 0.1;
y_CutOff= 3;

FontSize=8;
min_alpha = 0.01;
max_alpha = 0.5;
min_size = 1;
max_size = 75;


fh=figure('Name','Volcano Plot','Color','w','Tag','Volcano Plot figure','Units','inches');
fh.Position(3:4) = [2 2];
ah = axes(fh,'NextPlot','add','tag','Volcano Plot','box','off','Layer','top','FontSize',FontSize,'YAxisLocation','origin');
ah.LineWidth = 0.5;
switch PlotType
    case 'simple'
        ah.XGrid= 'off';
        ah.YGrid= 'off';
    otherwise
        ah.XGrid= 'on';
        ah.YGrid= 'on';
end



MaxVal_X_Significant = max(abs(x_data(y_data > y_CutOff)));

x_data(x_data>MaxVal_X_Significant) = MaxVal_X_Significant;
x_data(x_data<-MaxVal_X_Significant) = -MaxVal_X_Significant;

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

indx_pos_scatter = x_data_pos > x_CutOff & y_data_pos > y_CutOff;
indx_neg_scatter = x_data_neg < -x_CutOff & y_data_neg > y_CutOff;

cMap=colormap('bone');
cMap=flipud(cMap);
cMap = colormap(colorcet('L02','reverse',true));

DensScat(x_data_pos(~indx_pos_scatter),y_data_pos(~indx_pos_scatter),'TargetAxes',ah,'ColorBar',false,'ColorMap',cMap,'mSize',25);
sh_pos = scatter(ah,x_data_pos(indx_pos_scatter),y_data_pos(indx_pos_scatter),dist_pos_size(indx_pos_scatter),[1 0 0],'filled');
sh_pos.AlphaDataMapping = 'none';
sh_pos.AlphaData = dist_pos_alpha(indx_pos_scatter);
sh_pos.MarkerFaceAlpha = 'flat';

DensScat(x_data_neg(~indx_neg_scatter),y_data_neg(~indx_neg_scatter),'TargetAxes',ah,'ColorBar',false,'ColorMap',cMap,'mSize',25);
sh_neg = scatter(ah,x_data_neg(indx_neg_scatter),y_data_neg(indx_neg_scatter),dist_neg_size(indx_neg_scatter),[0 0 1],'filled');
sh_neg.AlphaDataMapping = 'none';
sh_neg.AlphaData = dist_neg_alpha(indx_neg_scatter);
sh_neg.MarkerFaceAlpha = 'flat';


xlabel(ah,xlabel_txt);
ylabel(ah,ylabel_txt);

min_x=min(x_data_neg,[],'all','omitnan');
max_x=max(x_data_pos,[],'all','omitnan');
nudge_x = (max_x-min_x)/20;

ah.XLim=[min_x-nudge_x max_x+nudge_x];

min_y=0;
max_y=max(y_data,[],'all','omitnan');
nudge_y = (max_y-min_y)/20;
ah.YLim=[0 max_y+nudge_y];

ah.YAxis.Label.HorizontalAlignment='center';
ah.YAxis.Label.VerticalAlignment='bottom';

switch PlotType
    case 'simple'
        pos_str=sprintf('n=%u',sum(indx_pos_scatter));
        text(ah,ah.XLim(2)-nudge_x,ah.YLim(2),pos_str,'Clipping','off','FontSize',FontSize,'HorizontalAlignment','right' );
        neg_str=sprintf('n=%u',sum(indx_neg_scatter));
        text(ah,ah.XLim(1)+nudge_x,ah.YLim(2),neg_str,'Clipping','off','FontSize',FontSize,'HorizontalAlignment','left' );


    otherwise
        for i = 2:length(ah.YAxis.TickValues)
            TickVal = ah.YAxis.TickValues(i);
            nPos = sum(y_data_pos>TickVal);
            pos_str=sprintf('n=%u',nPos);
            text(ah,ah.XLim(2)-nudge_x,TickVal,pos_str,'Clipping','off','FontSize',FontSize )
            nNeg = sum(y_data_neg>TickVal);
            neg_str=sprintf('n=%u',nNeg);
            text(ah,ah.XLim(1)+nudge_x,TickVal,neg_str,'Clipping','off','HorizontalAlignment','right','FontSize',FontSize)
        end
end
            ah.Position([1 3]) = [0.15 0.7];