function fh  = Methylation_Ctrl_Probe_BoxPlot(DATA,GroupVariableId, GroupsToUse, MatchedSamples)

RefDATA.IDs  = ["Restoration" "Staining.Green" "Staining.Red" "Extension.Green" "Extension.Red" "Hybridization.High.Medium" "Hybridization.Medium.Low" "Target.Removal.1" "Target.Removal.2" "Bisulfite.Conversion.I.Green" "Bisulfite.Conversion.I.Red" "Bisulfite.Conversion.II" "Specificity.I.Green" "Specificity.I.Red" "Specificity.II" "Non.polymorphic.Green" "Non.polymorphic.Red"];
RefDATA.Name = ["Restoration" "Staining green" "Staining red" "Extension green" "Extension red" "Hybridization high/medium" "Hybridization medium/low" "Target removal 1" "Target removal 2" "Bisulfite conversion I green" "Bisulfite conversion I red" "Bisulfite conversion II" "Specificity I green" "Specificity I red" "Specificity II" "Non-polymorphic green" "Non-polymorphic red"];
RefDATA.ThrVal = [0 5 5 5 5 1 1 1 1 1 1 1 1 1 1 5 5];


nVal = length(RefDATA.IDs);


fh = figure('Name','Control Probes boxplot','Color','w','Tag','Control Probes boxplot',...
    'Units','inches');
fh.Position(3:4) = [9 11];
th = tiledlayout(fh,3,6,'TileSpacing','compact','padding','compact');

for i=1:nVal
    ah = nexttile(th);
    [~, out] = PlotBoxPlotDATA(DATA,RefDATA.IDs{i},GroupVariableId,GroupsToUse,'MatchedSamples',MatchedSamples,...
        'CalcStats',true,'PlotStars',true,...
        'GroupColors',GetPalette('Tab10'),'MarkerSize',20,'MarkerLineWidth',0.5,'FontSize',8,'TargetAxes',ah);

    if out.Axes.YLim(1) < RefDATA.ThrVal(i)
        yline(ah,RefDATA.ThrVal(i),'color','r','Linewidth',1,'LineStyle','-');
    else
        text(ah,mean(out.Axes.XLim), out.Axes.YLim(1),sprintf('Threshold=%u',RefDATA.ThrVal(i)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8)
    end

end