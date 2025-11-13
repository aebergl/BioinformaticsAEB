function [fh RESULTS_DATA] = SampleDensityPanelFigure(DATA,IdColumn,Refsample,nRow,nCol,varargin)

FontSize = 10;
FigSize = [12 12];
TitleTxt = inputname(1);

ShowFigure = true;
FontSize = 10;
ExportPlot = false;
ExportDir = pwd;
MatchedSamplePairs = [];
ALLvsALL = false;
mSize = 10;
PointsToExclude = [];
if nargin > 5
    [varargin{:}] = convertStringsToChars(varargin{:});
end
ResolutionValue = 400;
LineWidth = 0.5;

DataType = "GeneExpression";
if DATA.nCol > 200000
    DataType = "Methylation";
end

% Check input
if nargin > 5
    ArgsList = {'FontSize','FigSize','TitleTxt','DisplayFigure','ExportPlot','ExportDir',...
        'MatchedSamplePairs','ALLvsALL','PointsToExclude'};
    for j=1:2:numel(varargin)
        ArgType = varargin{j};
        ArgVal = varargin{j+1};
        if ~strncmpi(ArgType,ArgsList,numel(ArgType))
            error('Invalid input option: %s', ArgType);
        else
            switch lower(ArgType)
                case 'fontsize'
                    FontSize = ArgVal;
                case 'figsize'
                    FigSize = ArgVal;
                case 'titletxt'
                    TitleTxt = ArgVal;
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
                case 'pointstoexclude'
                    PointsToExclude = ArgVal;
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
    if strcmpi('median',Refsample)
        x_ref =  median(DATA.X,1,"omitnan");
        OtherSampleIndx = 1:DATA.nRow;
        nArrays=DATA.nRow;
    else
        RefSampleIndx = strcmp(Refsample,SampleIds);
        OtherSampleIndx = find(~RefSampleIndx);
        nArrays=length(OtherSampleIndx);
        x_ref=DATA.X(RefSampleIndx,:);
    end
end

nPanes = nRow*nCol;
nImages = ceil(nArrays / nPanes);

fh = gobjects(nImages,1);
counter = 0;
switch DataType
    case 'GeneExpression'
        minX=min(DATA.X,[],'all')-0.3;
        maxX=max(DATA.X,[],'all')+0.3;
        VarNames = ["r Pearson" "r Spearman" "CCC" "n > 1.5-fold" "n > 2-fold" "% > 1.5-fold" "% > 2-fold"]';

    case 'Methylation'
        minX=min(DATA.X,[],'all')-0.01;
        maxX=max(DATA.X,[],'all')+0.01;
        VarNames = ["r Pearson" "r Spearman" "CCC" "n > 0.1" "n > 0.2" "% > 0.1" "% > 0.2"]';

end


nVar = length(VarNames);
RESULTS_DATA = CreateDataStructure(nArrays,nVar,[],[]);
% Add Info
RESULTS_DATA.Title = 'Sample Density';

RESULTS_DATA.ColId=VarNames;
RESULTS_DATA.RowAnnotation = strings(nArrays,1);
RESULTS_DATA.RowAnnotation(:,2) = TitleTxt;
RESULTS_DATA.RowAnnotationFields = {'Name','Title'};




for i = 1:nImages
    fh(i) = figure('Name','Array Image','Color','w','Tag','Array Image',...
        'Visible',ShowFigure,'Units','inches','Colormap',turbo);
    fh(i).Position(3:4) = FigSize;
    th = tiledlayout(fh(i),nRow,nCol,'TileSpacing','tight','padding','tight');
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
            DensScat(x_ref,y_sample, 'TargetAxes',ah,'AxisType','y=x','mSize',mSize,'PointsToExclude', PointsToExclude);
            xlabel(MatchedSamplePairs(counter,1),'FontSize',FontSize+2,'Interpreter','none');
            ylabel(MatchedSamplePairs(counter,2), 'FontSize',FontSize+2,'Interpreter','none')
            RESULTS_DATA.RowId(counter) = strcat(MatchedSamplePairs(counter,2)," vs. ",MatchedSamplePairs(counter,1));

        else
            y_sample=DATA.X(OtherSampleIndx(counter),:);
            DensScat(x_ref,y_sample, 'TargetAxes',ah,'AxisType','y=x','mSize',mSize,'PointsToExclude', PointsToExclude);
            xlabel(Refsample,'FontSize',FontSize+2,'Interpreter','none');
            ylabel(SampleIds(OtherSampleIndx(counter)), 'FontSize',FontSize+2,'Interpreter','none')
            RESULTS_DATA.RowId(counter) = SampleIds(OtherSampleIndx(counter));
            %RESULTS_DATA.RowId(counter) = strcat(MatchedSamplePairs(counter,2)," vs. ",Refsample);
        end
        nVal = sum(isreal(x_ref) & isreal(y_sample));
        ah.FontSize = FontSize;
        ah.XLim = [minX maxX];
        %ah.XTick = ceil(minX):2:floor(maxX);
        ah.YLim = [minX maxX];
        %ah.YTick = ceil(minX):2:floor(maxX);
        ah.XGrid ='on';
        ah.YGrid ='on';
        ah.Box='on';
        r_Pearson = corr(x_ref',y_sample','Type','Pearson','rows','pairwise');
        r_Spearman = corr(x_ref',y_sample','Type','Spearman','rows','pairwise');
        c = f_CCC([x_ref',y_sample'],0.05);
        switch DataType
            case 'GeneExpression'

                nDiff_1 = sum(abs(x_ref-y_sample) > 1);
                nDiff_2 = sum(abs(x_ref-y_sample) > 0.585);
                Str(1) = {[sprintf('rp=%.4f rs=%.4f ccc=%.4f',r_Pearson,r_Spearman,c{1}.est)]};
                Str(2) = {[sprintf('n > 2-fold: %u   n > 1.5-fold: %u',nDiff_1,nDiff_2)]};
                text(ah.XLim(1)+0.2,ah.YLim(2)+0.1,Str,'FontSize',FontSize,'VerticalAlignment','bottom','HorizontalAlignment','Left','Color','k');
                line([ah.XLim(1) ah.XLim(2)],[ah.YLim(1) ah.YLim(2)],'Linewidth',LineWidth,'LineStyle','-','color','k')
                line([ah.XLim(1) ah.XLim(2)-1],[ah.YLim(1)+1 ah.YLim(2)],'Linewidth',LineWidth,'LineStyle',':','color','k')
                line([ah.XLim(1)+1 ah.XLim(2)],[ah.YLim(1) ah.YLim(2)-1],'Linewidth',LineWidth,'LineStyle',':','color','k')
            case 'Methylation'
                nDiff_1 = sum(abs(x_ref-y_sample) > 0.1);
                nDiff_2 = sum(abs(x_ref-y_sample) > 0.2);
                Str(1) = {sprintf('rp=%.4f rs=%.4f ccc=%.4f',r_Pearson,r_Spearman,c{1}.est)};
                Str(2) = {sprintf('n > 0.1: %u   n > 0.2: %u',nDiff_1,nDiff_2)};
                text(ah.XLim(1)+0.01,ah.YLim(2)+0.01,Str,'FontSize',FontSize,'VerticalAlignment','bottom','HorizontalAlignment','Left','Color','k');
                line([ah.XLim(1) ah.XLim(2)],[ah.YLim(1) ah.YLim(2)],'Linewidth',LineWidth,'LineStyle','-','color','k')
                line([ah.XLim(1) ah.XLim(2)-0.1],[ah.YLim(1)+0.1 ah.YLim(2)],'Linewidth',LineWidth,'LineStyle',':','color','k')
                line([ah.XLim(1)+0.1 ah.XLim(2)],[ah.YLim(1) ah.YLim(2)-0.1],'Linewidth',LineWidth,'LineStyle',':','color','k')
                line([ah.XLim(1) ah.XLim(2)-0.2],[ah.YLim(1)+0.2 ah.YLim(2)],'Linewidth',LineWidth,'LineStyle','--','color','k')
                line([ah.XLim(1)+0.2 ah.XLim(2)],[ah.YLim(1) ah.YLim(2)-0.2],'Linewidth',LineWidth,'LineStyle','--','color','k')

        end

        RESULTS_DATA.X(counter,1) = r_Pearson;
        RESULTS_DATA.X(counter,2) = r_Spearman;
        RESULTS_DATA.X(counter,3) = c{1}.est;
        RESULTS_DATA.X(counter,4) = nDiff_1;
        RESULTS_DATA.X(counter,5) = nDiff_2;
        RESULTS_DATA.X(counter,4) = nDiff_1/nVal*100;
        RESULTS_DATA.X(counter,5) = nDiff_2/nVal*100;

        RESULTS_DATA.RowAnnotation(:,1) = RESULTS_DATA.RowId;
        drawnow
    end
    if TitleTxt
        title(th,TitleTxt,{'',''})
        th.Title.FontWeight = 'bold';
        th.Title.Interpreter = 'None';
        th.Title.FontSize = FontSize;
    end

    if ExportPlot
        if nImages > 1
            [~,name,ext] = fileparts(ExportPlot);
            tmpName = sprintf('%s_%u%s',name,i,ext);
            FullFileExport=fullfile(ExportDir,tmpName);
        else
            %FullFileExport=fullfile(ExportDir,ExportPlot);
            FullFileExport=fullfile(ExportDir,strcat(TitleTxt,".png"));
        end

        exportgraphics(th,FullFileExport,'Resolution',ResolutionValue)
        close(fh(i))
    end


end