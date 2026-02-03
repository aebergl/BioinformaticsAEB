
function [fh, pMV] = MethylationPlotMV(DATA,VariableIdentifier)



FontSize = 12;
FigSize = [5,5];
SamplId = false;
DataTipId = 'RowId';
SortOrder = 'chip';
SortOrder = 'pMV';
ChipBarcodeColumn = false;
ChipColumn = false;

% Get chip barcode
if ChipBarcodeColumn

else
    ChipBarcode = DATA.RowId;
end

% Get ChipId

if ChipColumn

else
    ChipId = extractBefore(ChipBarcode,'_');
end


% Calculate amount of missing values for each sample
nMV = sum(isnan(DATA.X),2);
pMV = nMV / DATA.nCol * 100;


switch lower(SortOrder)
    case 'chip'
        [~,SortIndx] = sort(ChipBarcode,'Ascend');
    case 'pmv'
        [~,SortIndx] = sort(ChipBarcode,'Descend');
    otherwise
        SortIndx = 1:DATA.nROW;
end


% Selection of X Y variable id
if VariableIdentifier
    indx = strcmpi(VariableIdentifier,DATA.ColAnnotationFields);
    if ~any(indx)
        error('Error. \n%s not found in DATA.ColAnnotationFields',VariableIdentifier);
    elseif sum(indx) > 1
        error('Warning. \nMultiple matches for %s found in DATA.ColAnnotationFields',VariableIdentifier);
    else
        DATA_ID = DATA.ColAnnotation(:,indx);
    end
else
    DATA_ID = DATA.ColId;
end


switch lower(DataTipId)
    case 'rowid'
        SampleId = DATA.RowId;
    otherwise
        indx_SampleId = strcmpi(DataTipId,DATA.RowAnnotationFields);
        if ~any(indx_SampleId)
            error('Error. \n%s not found in DATA.RowAnnotationFields',DataTipId);
        elseif sum(indx_VarId) > 1
            error('Warning. \nMultiple matches for %s found in DATA.RowAnnotationFields',DataTipId);
        else
            SampleId = DATA.RowAnnotation(:,indx_SampleId);
        end
end

% Sort Data

nMV = nMV(SortIndx);
pMV = pMVorig(SortIndx);




% Percent Missing values

fh=figure('Name','Bar Plot','Color','w','Tag','Bar Plot','Units','inches');
fh.Position(3:4) = FigSize;
ah = axes(fh,'NextPlot','add','tag','Scatter Plot','box','on','Layer','top','FontSize',FontSize,'YGrid','on');
ah.LineWidth = 0.5;
bh = bar(ah,pMV,'FaceColor','flat');
ylabel('% Missing Values');
ah.XTick = 1:DATA.nRow;
ah.XLim = [0.5 DATA.nRow + 0.5];
ah.XTickLabel = ChipId;


            row = dataTipTextRow('',DATA_n200_cg.RowId);
            bh.DataTipTemplate.DataTipRows = row;
            bh.DataTipTemplate.Interpreter='none';

bh.FaceColor='flat'
indx_orange=bh.YData>5;
indx_red=bh.YData>10;
bh.CData(indx_orange,:) = repmat(GetPalette('Tab10',2),sum(indx_orange),1);
bh.CData(indx_red,:) = repmat(GetPalette('Tab10',4),sum(indx_red),1);
yline(ah,[5 5],'color','k','Linewidth',1)
yline(ah,[10 10],'color','k','Linewidth',1)
