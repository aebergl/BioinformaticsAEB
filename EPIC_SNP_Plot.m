function [fh] = EPIC_SNP_Plot(DATA,SortingVariable,SimilarityThreshold,FigSize)

SortingDirection = 'ascend';

FontSize = 7;
minSize = 20;
maxSize = 100;
LineWidth = 0.5;
GridLines = 'on';
RightMargin = 0.5;




if isempty(SortingVariable)
    SortIndx = 1:DATA.nRow;
    indx_SortingVar = 1;
else
    indx_SortingVar = strcmpi(SortingVariable,DATA.RowAnnotationFields);
    if any(indx_SortingVar)
        [~,SortIndx] = sort(DATA.RowAnnotation(:,indx_SortingVar));
    else
        error('%s not found',ColorType)
    end
end

X_data = DATA.X(SortIndx,:);
SampleId = DATA.RowId(SortIndx);
if isempty(SortingVariable)
    SampleAnnotation = SampleId;
else
    SampleAnnotation = DATA.RowAnnotation(SortIndx,:);
end


if size(X_data,2) < 100
    if anynan(X_data)
        x = pdist(X_data,@naneucdist);
        
    else
x = pdist(X_data,'euclidean');
    end
    DataType = "SNP";
else
    x=zeros((DATA.nRow*(DATA.nRow-1)/2),1);
    counter = 0;
    for i=1:DATA.nRow-1
        for j=i+1:DATA.nRow
            counter = counter + 1;
            x(counter) = sum(abs(X_data(i,:)-X_data(j,:)) > 0.5);
        end
    end
    DataType = "Methylation";
    x=log10(x+1);
end

fh_H=figure('Name','Distance Histogram','Color','w','Tag','Distance Histogram','Units','inches');
fh_H.Position(3:4) = FigSize;
ah_H = axes(fh_H,'NextPlot','add','tag','Volcano Plot','box','on','Layer','top','FontSize',FontSize,...
    'PositionConstraint','outerposition','Clipping','off');
ah_H.LineWidth = LineWidth;
histogram(ah_H,x,100)
xlabel('Similarity Score');
ylabel('Count');


fh=figure('Name','Distance Histogram','Color','w','Tag','Distance Histogram','Units','inches');
fh.Position(3:4) = FigSize;
ah = axes(fh,'NextPlot','add','tag','Volcano Plot','box','on','Layer','top','FontSize',FontSize,...
    'PositionConstraint','outerposition','Clipping','off','YDir','Reverse');
ah.LineWidth = LineWidth;


x = squareform(x);
fh=imagesc(ah,x);

axis equal
fh.Parent.YLim =[0.5 DATA.nRow+0.5];
fh.Parent.XLim =[0.5 DATA.nRow+0.5];
fh.Parent.XTick=[];
fh.Parent.YTick=[];
fh.Parent.YTick=1:DATA.nRow;
fh.Parent.YTickLabel = SampleAnnotation(:,indx_SortingVar);
fh.Parent.YAxis.TickLabelInterpreter='none';
colorbar




switch DataType
    case "SNP"
        Xmask=ones(size(x))*inf;

        Xsim = x+triu(Xmask);
        [ia, ib] = find(Xsim < SimilarityThreshold);
    case "Methylation"
        Xmask=ones(size(x))*inf;

        Xsim = x+triu(Xmask);
        [ia, ib] = find(Xsim < SimilarityThreshold);
end

numSimilar = length(ia);

[numId,numMatching] = GroupCount([ia;ib],0);

fprintf("Barcode A\tSample Id A\tnMatching A\tBarcode B\tSample Id B\tnMatching B\tSimilarity Score\n")
for i=1:numSimilar

    numA=numMatching(ia(i)==numId);
    numB=numMatching(ib(i)==numId);
    fprintf('%s\t%s\t%u\t%s\t%s\t%u\t%f\n',SampleId{ia(i)},SampleAnnotation{ia(i),indx_SortingVar},numA,SampleId{ib(i)},SampleAnnotation{ib(i),indx_SortingVar},numB,x(ia(i),ib(i)))

end

end

function D2 = naneucdist(XI,XJ)  
%NANEUCDIST Euclidean distance ignoring coordinates with NaNs
n = size(XI,2);
sqdx = (XI-XJ).^2;
nstar = sum(~isnan(sqdx),2); % Number of pairs that do not contain NaNs
nstar(nstar == 0) = NaN; % To return NaN if all pairs include NaNs
D2squared = sum(sqdx,2,'omitnan').*n./nstar; % Correction for missing coordinates
D2 = sqrt(D2squared);

end
