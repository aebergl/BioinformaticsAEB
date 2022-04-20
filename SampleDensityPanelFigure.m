function fh = SampleDensityPanelFigure(DATA,IdColumn,Refsample,nRow,nCol,varargin)


PageWidth = 7;
ShowFigure = true;
FontSize = 10;
ExportPlot = false;
ExportDir = pwd;

if nargin > 5
    [varargin{:}] = convertStringsToChars(varargin{:});
end


% Check input
if nargin > 5
    ArgsList = {'DisplayFigure', 'ExportPlot','ExportDir'};
    for j=1:2:numel(varargin)
        ArgType = varargin{j};
        ArgVal = varargin{j+1};
        if ~strncmpi(ArgType,ArgsList,numel(ArgType))
            error('Invalid input option: %s', ArgType);
        else
            switch lower(ArgType)
                case 'displayfigure'
                    ShowFigure = ArgVal;
                case 'exportplot'
                    ExportPlot = ArgVal;
                case 'exportdir'
                    ExportDir = ArgVal;

            end
        end
    end
end


if isempty(IdColumn)
    SampleIds = DATA.RowId
else
    IDColumnIndx = strcmp(IdColumn,DATA.RowAnnotationFields);
    SampleIds = DATA.RowAnnotation(:,IDColumnIndx);
end

RefSampleIndx = strcmp(Refsample,SampleIds);

OtherSampleIndx = find(~RefSampleIndx);

nArrays=length(OtherSampleIndx);
nPanes = nRow*nCol;
nImages = ceil(nArrays / nPanes);

fh = gobjects(nImages,1);
counter = 0;
x_ref=DATA.X(RefSampleIndx,:);
minX=min(DATA.X,[],'all')-0.1;
maxX=max(DATA.X,[],'all')+0.1;

for i = 1:nImages
    fh(i) = figure('Name','Array Image','Color','w','Tag','Array Image',...
        'Visible',ShowFigure,'Units','inches','Colormap',turbo);
    fh(i).Position(3) = PageWidth;
    fh(i).Position(4) = PageWidth*nRow/nCol;
    th = tiledlayout(fh(i),nRow,nCol,'TileSpacing','tight','padding','compact');
    %th.Title.String = Refsample;
    %th.Title.FontWeight = 'bold';
    %th.Title.Interpreter = 'None';
    %th.Title.FontSize = FontSize;
    CurrentPane = 0;
    while CurrentPane < nPanes && counter < nArrays
        counter = counter + 1;
        CurrentPane = CurrentPane + 1;
        ah = nexttile(th);
        y_sample=DATA.X(OtherSampleIndx(counter),:);
        DensScat(x_ref,y_sample, 'TargetAxes',ah,'AxisType','y=x','mSize',5);

        xlabel(Refsample,'FontSize',FontSize+2,'Interpreter','none');
        ylabel(SampleIds(OtherSampleIndx(counter)), 'FontSize',FontSize+2,'Interpreter','none')
        ah.XLim = [minX maxX];
        ah.XTick = ceil(minX):2:floor(maxX);
        ah.YLim = [minX maxX];
        ah.YTick = ceil(minX):2:floor(maxX);
        ah.XGrid ='on';
        ah.YGrid ='on';
        ah.Box='on';
        r_corr = corr(x_ref',y_sample','rows','pairwise');
        %r_corr=1.0;
        nDiff_2 = sum(abs(x_ref-y_sample) > 1);
        nDiff_15 = sum(abs(x_ref-y_sample) > 0.585);
        %Str(1) = {['\itn\rmFC>2=',sprintf('%u',nDiff_2)]};
        %Str(1) = {['\itn\rmFC>2=',sprintf('%u',nDiff_2)]};
        Str(1) = {['r = ',sprintf('%.5f',r_corr)]};
        Str(2) = {['Num > 2-fold: ',sprintf('%u',nDiff_2)]};
        Str(3) = {['Num > 1.5-fold: ',sprintf('%u',nDiff_15)]};
        

        text(ah.XLim(1)+0.2,ah.YLim(2)+0.1,Str,'FontSize',FontSize,'VerticalAlignment','bottom','HorizontalAlignment','Left','Color','k');

        line([ah.XLim(1) ah.XLim(2)-1],[ah.YLim(1)+1 ah.YLim(2)],'Linewidth',1,'LineStyle','-.','color','r')
        line([ah.XLim(1)+1 ah.XLim(2)],[ah.YLim(1) ah.YLim(2)-1],'Linewidth',1,'LineStyle','-.','color','r')


    end
    if ExportPlot
        if nImages > 1
            [~,name,ext] = fileparts(ExportPlot);
            tmpName = sprintf('%s_%u%s',name,i,ext);
            FullFileExport=fullfile(ExportDir,tmpName);
        else
            FullFileExport=fullfile(ExportDir,ExportPlot);
        end

        exportgraphics(th,FullFileExport,'Resolution',200)
    end


end