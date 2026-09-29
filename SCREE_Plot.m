function SCREE_Plot(PCAModel,TitleText,Type)

figure
numValues = PCAModel.NumComp;
if Type == 1
    plot(PCAModel.ExplVar,'o','Linewidth',2,'MarkerSize',10,'Color','k')
    hold on
    plot(PCAModel.ExplVar,'-','Linewidth',2,'MarkerSize',10,'Color','r')
    h=gca;
    h.LineWidth=1;
    set(gca,'FontSize',16);
    set(gca,'Xtick',1:1:numValues)
    xlabel('PCA component','FontSize',18,'Interpreter','none')
    ylabel('Explained variation (%)','FontSize',18,'Interpreter','none')
elseif Type == 2
        plot([0; PCAModel.ExplVarCum],'o','Linewidth',2,'MarkerSize',10,'Color','k')
        hold on
        plot([0; PCAModel.ExplVarCum],'-','Linewidth',2,'MarkerSize',10,'Color','r')
        h=gca;
        h.LineWidth=1;
        set(gca,'FontSize',16);
        set(gca,'Xtick',0:10:numValues)
        set(gca,'XLim',[0 numValues])
        xlabel('PCA component','FontSize',18,'Interpreter','none')
        ylabel('Explained variation (%)','FontSize',18,'Interpreter','none')
elseif Type == 3
        plot(PCAModel.Eig,'o','Linewidth',2,'MarkerSize',10,'Color','k')
        hold on
        plot(PCAModel.Eig,'-','Linewidth',2,'MarkerSize',10,'Color','r')
        h=gca;
        h.LineWidth=1;
        set(gca,'FontSize',16);
        set(gca,'Xtick',1:1:numValues)
        xlabel('PCA component','FontSize',18,'Interpreter','none')
        ylabel('Eigenvalue','FontSize',18,'Interpreter','none')

end
    
    
    title(TitleText,'FontSize',18,'Interpreter','none')
    set(gcf, 'Color', 'w');