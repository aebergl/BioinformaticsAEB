function h = MethylationProbeCorrPlot(DATA,probes,varargin);


FontSize = 7;
%CorrType='Pearson';
CorrType='Spearman';
FigureSize = [5 6];
AxesSize = [4.2 4.2];
LeftBorderWidth = 1.6;
BottomeBorderHight = 0.8;
[CMap, descriptorname, description] = colorcet('D1');


[~,ia,ib]=intersect(probes,DATA.ColId,'Stable');
probes = probes(ia);

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
h.InnerPosition = [LeftBorderWidth BottomeBorderHight AxesSize];
h.FontSize = FontSize;
h.YDisplayLabels = probes;

h.XDisplayLabels = strings(N,1);

ah = axes(h.Parent,'Units','inches','NextPlot','add','Box','on',...
    'FontSize',FontSize,'Linewidth',0.5,'XGrid','on','YGrid','on','PositionConstraint','Outerposition',...
    'InnerPosition',[LeftBorderWidth  h.InnerPosition(2)+AxesSize(1)+0.15    AxesSize]);
scatter(ah,repmat(1:N,M,1),X,'XJitter','density','XJitterWidth',0.5 )
ah.XLim=[0.5 N+0.5];
ah.XTick = 1:N;
ah.XAxisLocation='bottom';
ah.XTickLabelRotation = 0;
ylabel(ah,'\beta-value');


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
%%x_str = sprintf('\\fontsize{16} \\bf Gene: {\\color[rgb]%sTSS1500 \\color[rgb]%sTSS200 \\color[rgb]%s5UTR \\color[rgb]%s1stExon \\color[rgb]%sBody \\color[rgb]%s3UTR}',S_C1,S_C2,S_C3,S_C4,S_C5,S_C6);
x_str = sprintf('\\fontsize{14}\\bf{\\color[rgb]%sTSS1500\n\\color[rgb]%sTSS200\n\\color[rgb]%s5UTR\n\\color[rgb]%s1stExon\n\\color[rgb]%sBody\n\\color[rgb]%s3UTR}',S_C1,S_C2,S_C3,S_C4,S_C5,S_C6);
x_str = sprintf('\\fontsize{8}\\bf{\\color[rgb]%sTSS1500\n\\color[rgb]%sTSS200\n\\color[rgb]%s5UTR\n\\color[rgb]%s1stExon\n\\color[rgb]%sExonBnd\n\\color[rgb]%s3UTR}',S_C1,S_C2,S_C3,S_C4,S_C5,S_C6);

Ypos = 0;
% GeneAccession_indx = ~isundefined(S.GeneAccession);
% for j = 1:S.numTranscripts
%     Ypos = Ypos - 0.05;
%    for i=1:numProbes
%        GeneBodyType = cellstr(S.GeneGroup(i,j));
%        indx = strcmp(GeneBodyType,GeneGroupTypes);
%        if any(indx)
%            line([ i-0.5 i+0.5], [Ypos Ypos],'LineStyle','-','LineWidth',4,'Color',GeneColor(indx,:))
%        end
%        indx = find(GeneAccession_indx(:,j));
%        text(0.5,Ypos,cellstr(S.GeneAccession(indx(1),j)),'FontSize',6,'HorizontalAlignment','right','Interpreter','none')
% 
%    end
% end
% line([ 0.5 numProbes + 0.5], [0 0],'LineStyle','-','LineWidth',1,'Color',[0 0 0]);







