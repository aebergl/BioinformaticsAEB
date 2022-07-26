function fh = SampleDensityPanelFigure(DATA,IdColumn,Refsample,nRow,nCol,varargin)


PageWidth = 7;
ShowFigure = true;
FontSize = 10;
ExportPlot = false;
ExportDir = pwd;
MatchedSamplePairs = [];
ALLvsALL = false;
if nargin > 5
    [varargin{:}] = convertStringsToChars(varargin{:});
end


% Check input
if nargin > 5
    ArgsList = {'DisplayFigure', 'ExportPlot','ExportDir','MatchedSamplePairs','ALLvsALL'};
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
                case 'matchedsamplepairs'
                    MatchedSamplePairs = ArgVal;
                case 'allvsall'
                    ALLvsALL = ArgVal;

            end
        end
    end
end

if isempty(IdColumn)
    SampleIds = DATA.RowId;
else
    IDColumnIndx = strcmp(IdColumn,DATA.RowAnnotationFields);
    SampleIds = DATA.RowAnnotation(:,IDColumnIndx);
end
if ALLvsALL
    n=length(SampleIds);
    MatchedSamplePairs = cell(n*(n-1)/2,2);
    counter = 0;
    for i = 1:n-1
        for j = i+1:n
            counter = counter + 1;
            MatchedSamplePairs(counter,1) = SampleIds(i);
            MatchedSamplePairs(counter,2) = SampleIds(j);
        end
    end

end

if ~isempty(MatchedSamplePairs)
    nArrays=size(MatchedSamplePairs,1);
else
    RefSampleIndx = strcmp(Refsample,SampleIds);
    OtherSampleIndx = find(~RefSampleIndx);
    nArrays=length(OtherSampleIndx);
    x_ref=DATA.X(RefSampleIndx,:);
end
nPanes = nRow*nCol;
nImages = ceil(nArrays / nPanes);

fh = gobjects(nImages,1);
counter = 0;

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
        if ~isempty(MatchedSamplePairs)
            RefSampleIndx_x = strcmp(MatchedSamplePairs(counter,1),SampleIds);
            RefSampleIndx_y = strcmp(MatchedSamplePairs(counter,2),SampleIds);
            x_ref=DATA.X(RefSampleIndx_x,:);
            y_sample=DATA.X(RefSampleIndx_y,:);
            DensScat(x_ref,y_sample, 'TargetAxes',ah,'AxisType','y=x','mSize',5);
            xlabel(MatchedSamplePairs(counter,1),'FontSize',FontSize+2,'Interpreter','none');
            ylabel(MatchedSamplePairs(counter,2), 'FontSize',FontSize+2,'Interpreter','none')

        else
            y_sample=DATA.X(OtherSampleIndx(counter),:);
            DensScat(x_ref,y_sample, 'TargetAxes',ah,'AxisType','y=x','mSize',5);
            xlabel(Refsample,'FontSize',FontSize+2,'Interpreter','none');
            ylabel(SampleIds(OtherSampleIndx(counter)), 'FontSize',FontSize+2,'Interpreter','none')
        end
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