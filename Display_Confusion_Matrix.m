function [fh] = Display_Matrix(C,FigSize,TickLabels,XLabelTxt, YLabelTxt)

p.LineWidth = 0.5;
p.LineStyle = '-';
p.MatrixFontSize = 10;
p.FontSize = 10;

fh=figure('Name','HeatMap','Color','w','Tag','HeatMap figure','Units','inches');
fh.Position(3:4) = FigSize;
ah = axes(fh,'NextPlot','add','tag','HeatMap','Units','inches',...
    'Box','On','LineWidth',p.LineWidth,'ydir', 'reverse');
ah.FontSize = p.FontSize;
ah.PositionConstraint = "innerposition";

[nRows,nCols] = size(C);
[px,py] = GetPatch(nRows,nCols);

ph = patch(ah,px, py, reshape(C,1,nRows*nCols), 'linestyle', p.LineStyle);
colormap(ah,'sky');
cmap = colormap;

ah.XLim = [.5 nCols+.5];
ah.YLim = [.5 nRows+.5];

Cp = C ./ sum(C,'all') * 100;

hold on;
for r = 1:nRows
    for c = 1:nCols
        txt{1} = num2str(C(r,c));
        txt{2} = sprintf('(%.1f%%)',Cp(r,c));
        text(ah,c, r, txt, ...
            "FontSize",p.MatrixFontSize, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'Color', 'black', 'FontWeight', 'bold');
    end
end
hold off;
xlabel(XLabelTxt);
ylabel(YLabelTxt);
ah.XTick = 1:nCols;
ah.YTick = 1:nRows;

ah.XTickLabel = TickLabels;
ah.YTickLabel = TickLabels;

ah.Position(1:2) = [0.8 0.5];
ah.Position(3:4) = [FigSize(1)-ah.Position(1)-0.05, FigSize(1)-ah.Position(1)-0.05];


function [px,py] = GetPatch(nRows,nCols)
%Calculate Patch coordinates
%https://stackoverflow.com/questions/6614207/how-to-export-non-blurry-eps-images
px = bsxfun(@plus, [-0.5; 0.5; 0.5; -0.5], reshape(1:nCols, [1 1 nCols]));
py = bsxfun(@plus, [-0.5; -0.5; 0.5; 0.5], 1:nRows);

px = reshape(repmat(px, [1 nRows 1]), 4, nCols*nRows);
py = reshape(repmat(py, [1 1 nCols]), 4, nCols*nRows);

end


end
