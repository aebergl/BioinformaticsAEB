function fh = DensityScatterMatrix(DATA)

PageWidth = 10;
FontSize=8;
Intepreter = 'none';

X = DATA.X;
[N,K] = size(X);

fh = figure('Name','Array Image','Color','w','Tag','Array Image',...
    'Units','inches');
fh.Position(3) = PageWidth;
fh.Position(4) = PageWidth;


th = tiledlayout(fh,N,N,'TileSpacing','none','padding','tight','TileIndexing','columnmajor','PositionConstraint','outerposition');

% A = true(N,N);
% indx_lower = find(tril(A,-1));
% counter = 0;
for i=1:N-1
    for j=i+1:N
        nTile = sub2ind([ N N],j,i);
        ah = nexttile(nTile);
        DensScat(X(i,:),X(j,:),'TargetAxes',ah,'AxisType','y=x','ColorBar',false,'mSize',10);
        ah.XTick=[];
        ah.YTick=[];
        ah.Box='on';
        if i==1
            ylabel(DATA.RowId(j),'FontSize',FontSize,'Interpreter',Intepreter);
        end
        if j==N
            xlabel(DATA.RowId(i),'FontSize',FontSize,'Interpreter',Intepreter);
        end

        nTile = sub2ind([ N N],i,j);
         ah = nexttile(nTile);
         axis equal;
        ah.XLim=[-1 1];
        ah.YLim=[-1 1];
        ah.XTick=[];
        ah.YTick=[];
        ah.Box='on';

        [rs,ps] =  corr(X(i,:)',X(j,:)','rows','pairwise','type','Spearman');
        text(ah,0,0,sprintf('%.3g',rs),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',FontSize+4,'Interpreter',Intepreter)
        drawnow
    end
end
for i=1:N
    nTile = sub2ind([ N N],i,i);
    ah = nexttile(nTile);
    axis equal;
    ah.XLim=[-1 1];
    ah.YLim=[-1 1];
    ah.XTick=[];
    ah.YTick=[];
    ah.Box='on';
    text(ah,0,0,DATA.RowId{i},'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',FontSize+2,'Interpreter',Intepreter)

end