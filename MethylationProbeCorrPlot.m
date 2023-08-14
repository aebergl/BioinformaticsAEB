function h = MethylationProbeCorrPlot(DATA,probes,varargin);


FontSize = 8;
CorrType='Pearson';
%CorrType='Spearman';
FigureSize = [5 6];
AxesSize = [2 2];
LeftBorderWidth = 1.6;
BottomeBorderHight = 0.8;
[CMap, descriptorname, description] = colorcet('D1');


% 
% if length(probes) > 1
    [~,ia,ib]=intersect(probes,DATA.ColId,'Stable');
    probes = probes(ia);
% else
%     GeneSingleColumn = find(strcmp('UCSC_REFGENE_NAME_SINGLE',DATA.ProbeAnnotationColumns));
%     GeneIdColumn = find(strcmp('UCSC_REFGENE_NAME',DATA.ProbeAnnotationColumns));
%     ChrPosColumn = find(strcmp('MAPINFO',DATA.ProbeAnnotationColumns));
% 
%     WindowIds = DATA.ProbeAnnotation(:,find(strcmp('UCSC_REFGENE_NAME_SINGLE',DATA.ProbeAnnotationColumns)));
%     WindowIdsAll = DATA.ProbeAnnotation(:,find(strcmp('UCSC_REFGENE_NAME',DATA.ProbeAnnotationColumns)));
% 
%     ProbeIndx = find(strcmp(probes,WindowIds));
%     ProbeIndxMultiple_1 = find(contains(WindowIdsAll,[probes,';']));
%     ProbeIndxMultiple_2 = find(contains(WindowIdsAll,[';',probes]));
%     if ~isempty(ProbeIndxMultiple_1)
%         ProbeIndx = union(ProbeIndx,ProbeIndxMultiple_1);
%     end
%     if ~isempty(ProbeIndxMultiple_2)
%         ProbeIndx = union(ProbeIndx,ProbeIndxMultiple_2);
%     end
% 
%     numProbes = length(ProbeIndx);
%     % Check for righ gene
%     rem_indx=[];
%     for i=1:numProbes
%         Gene=DATA.ProbeAnnotation{ProbeIndx(i),GeneIdColumn};
%         tmp = strsplit(Gene,';');
%         if sum(strcmp(probes,tmp)) == 0
%             rem_indx = [rem_indx;i];
%         end
%     end
%     ProbeIndx(rem_indx)=[];
%     Chr_pos_text = DATA.ProbeAnnotation(ProbeIndx,ChrPosColumn);
%     Chr_pos = sscanf(sprintf('%s,',Chr_pos_text{:}),'%i,');
%     [Chr_pos,sort_indx] = sort(Chr_pos,'ascend');
%     ProbeIndx = ProbeIndx(sort_indx);
%     GenId = probes;
%     probes = DATA.ProbeId(ProbeIndx);
% end

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
% ah = axes(fh,'NextPlot','add','tag','Gene Sample Plot','Box','on','FontSize',FontSize,'Linewidth',0.5,...
%     'ActivePositionProperty','outerposition','XGrid','on','YGrid','on');
% ah.LineWidth = 0.5;
% ah.Colormap=CMap;






