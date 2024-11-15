function h = MethylationGeneModelPlot(DATA,GeneModel,GeneName,varargin)

FontSize = 10;
%CorrType='Pearson';
CorrType='Spearman';
FigureSize = [5 6];
AxesSize = [4 4];
LeftBorderWidth = 1.6;
BottomBorderHight = 1;

[CMap, ~, ~] = colorcet('D1');

GENE = GeneModel.GeneData(strcmp(GeneName,GeneModel.GeneList));
if isempty(GENE)
    error('%s not found in the data!',GeneName)
end


[~,ia,ib]=intersect(GENE.probes,DATA.ColId,'Stable');
probes = GENE.probes(ia);

X = DATA.X(:,ib);
    if strcmpi(CorrType,'Pearson')
        CorrMat = corr(X,'type','Pearson','rows','pairwise');
    elseif strcmpi(CorrType,'Spearman')
        CorrMat = corr(X,'type','Spearman','rows','pairwise');
    end
% if flipMat
%     CorrMat = flipud(CorrMat);
%     CorrMat = fliplr(CorrMat);
%     probes = flip(probes);
% end

[M,N] = size(X);

fh = figure('Name','Probe correlation heatmap','Color','w','Tag','Probe correlation heatmap','Units','inches');
fh.Position(3:4) = FigureSize;

h=heatmap(fh,CorrMat,'ColorLimits',[-1 1]);
h.Colormap = CMap;
h.Units='inches';
h.PositionConstraint = 'innerposition';
h.InnerPosition = [LeftBorderWidth BottomBorderHight AxesSize];
h.FontSize = FontSize;
h.YDisplayLabels = probes;

h.XDisplayLabels = strings(N,1);


% Scatter Plot
ah_SC = axes(h.Parent,'Units','inches','NextPlot','add','Box','on',...
    'FontSize',FontSize,'Linewidth',0.5,'XGrid','on','YGrid','on','PositionConstraint','Outerposition',...
    'InnerPosition',[LeftBorderWidth  h.InnerPosition(2)+AxesSize(1)+0.25    AxesSize]);
% scatter(ah_SC,repmat(1:N,M,1),X,'XJitter','density','XJitterWidth',0.5 )
% ah_SC.XLim=[0.5 N+0.5];
% ah_SC.XTick = 1:N;
% ah_SC.XAxisLocation='bottom';
% ah_SC.XTickLabelRotation = 0;
Group=DATA.RowAnnotation(:,4);
Group = repmat(Group,N,1);
X_pos = repmat(1:N,M,1);
boxchart(ah_SC,X_pos(:),X(:),'GroupByColor',Group)
ah_SC.XLim=[0.5 N+0.5];
ah_SC.XTick = 1:N;
ah_SC.XAxisLocation='bottom';
ah_SC.XTickLabelRotation = 0;

ylabel(ah_SC,'\beta-value');
title(GeneName,'FontWeight','normal')

% Gene Group Plot

GeneGroupTypes = {'TSS1500','TSS200','5''UTR','1stExon','ExonBnd','3''UTR'};
%GeneGroupTypes = {'tss_1500','tss_200','tss_body'};
%GeneGroupTypes = {'TSS1500','TSS200','5''UTR','1stExon','Body','3''UTR'};

GeneColor = GetPalette('aeb01',[3 6 7 8 9 11]);
S_C1 = sprintf('{%f, %f, %f}',GeneColor(1,1),GeneColor(1,2),GeneColor(1,3));
S_C2 = sprintf('{%f, %f, %f}',GeneColor(2,1),GeneColor(2,2),GeneColor(2,3));
S_C3 = sprintf('{%f, %f, %f}',GeneColor(3,1),GeneColor(3,2),GeneColor(3,3));
S_C4 = sprintf('{%f, %f, %f}',GeneColor(4,1),GeneColor(4,2),GeneColor(4,3));
S_C5 = sprintf('{%f, %f, %f}',GeneColor(5,1),GeneColor(5,2),GeneColor(5,3));
S_C6 = sprintf('{%f, %f, %f}',GeneColor(6,1),GeneColor(6,2),GeneColor(6,3));
%x_str = sprintf('\\fontsize{16} \\bf Gene: {\\color[rgb]%sTSS1500 \\color[rgb]%sTSS200 \\color[rgb]%s5UTR \\color[rgb]%s1stExon \\color[rgb]%sBody \\color[rgb]%s3UTR}',S_C1,S_C2,S_C3,S_C4,S_C5,S_C6);
%x_str = sprintf('\\fontsize{14}\\bf{\\color[rgb]%sTSS1500\n\\color[rgb]%sTSS200\n\\color[rgb]%s5UTR\n\\color[rgb]%s1stExon\n\\color[rgb]%sBody\n\\color[rgb]%s3UTR}',S_C1,S_C2,S_C3,S_C4,S_C5,S_C6);
x_str = sprintf('\\fontsize{8}\\bf{\\color[rgb]%sTSS1500\n\\color[rgb]%sTSS200\n\\color[rgb]%s5UTR\n\\color[rgb]%s1stExon\n\\color[rgb]%sExonBnd\n\\color[rgb]%s3UTR}',S_C1,S_C2,S_C3,S_C4,S_C5,S_C6);

ah_GG = axes(h.Parent,'Units','inches','NextPlot','add','Box','on','Clipping','off',...
    'FontSize',FontSize,'Linewidth',0.5,'XGrid','off','YGrid','off','PositionConstraint','Outerposition',...
    'InnerPosition',[LeftBorderWidth  0.05  AxesSize(1) BottomBorderHight-0.1]);

GeneAccession_indx = ~isundefined(GENE.GeneAccession);
Ypos = 0;
ah_GG.XLim=[0.5 N+0.5];
for j = 1:GENE.numTranscripts
    Ypos = Ypos + 1;
   for i=1:length(probes)
       indx_probe = strcmpi(probes{i},GENE.probes);
       GeneBodyType = cellstr(GENE.GeneGroup(indx_probe,j));
       indx = strcmp(GeneBodyType,GeneGroupTypes);
       if any(indx)
           line([ i-0.5 i+0.5], [Ypos Ypos],'LineStyle','-','LineWidth',4,'Color',GeneColor(indx,:))
       end
       indx = find(GeneAccession_indx(:,j));
       text(-0.1,Ypos,cellstr(GENE.GeneAccession(indx(1),j)),'FontSize',FontSize,'HorizontalAlignment','right','Interpreter','none')

   end
end
%line([ 0.5 N + 0.5], [0 0],'LineStyle','-','LineWidth',1,'Color',[0 0 0]);
text(ah_GG,N+1,GENE.numTranscripts+0.5,x_str,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',FontSize)

 
 ah_GG.XTick = [];
 ah_GG.YTick = [];
 ah_GG.YLim=[0 GENE.numTranscripts+0.5];







