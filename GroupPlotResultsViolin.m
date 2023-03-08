function fh = GroupPlotResultsViolin(DATA,GroupId,Groups,Y_Type,Y_Val_CutOff,p_Val_Type,p_Val_Cutoff,varargin)
FontSize = 10;
Plot_ALL=false;
TopNValues = 0;
TopPrctile = 0;

i=0;

while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'TopNValues')
        i = i + 1;
        TopNValues = varargin{i};
    elseif strcmpi(varargin{i},'TopPrctile')
        i = i + 1;
        TopPrctile = varargin{i};


    end
end


% if ~iscell(Groups)
%     Groups = {Groups}
% end


nVals = DATA.nRow;

GroupColumn = strcmp(GroupId,DATA.RowAnnotationFields);
if ~any(GroupColumn)
    error('No group Column Found')
end
nGroups = numel(Groups);



numValuesPerGroup = zeros(nGroups,1);
indx_Groups=false(DATA.nRow,nGroups);
GroupVariable = cell(DATA.nRow,1);
GroupName = cell(nGroups,1);

for i=1:nGroups
    tmp_name = Groups{i};
    if ~iscell(tmp_name)
        tmp_name = {tmp_name};
    end
    SampleIndxMatrix = false(DATA.nRow,numel(tmp_name));
    for j=1:length(tmp_name)
        SampleIndxMatrix(:,j) = strcmp(tmp_name{j},DATA.RowAnnotation(:,GroupColumn));
    end
    indx_Groups(:,i) = any(SampleIndxMatrix,2);
    GroupVariable(indx_Groups(:,i)) = tmp_name(1);
    numValuesPerGroup(i) = sum(indx_Groups(:,i));
    GroupName(i) = tmp_name(1);
end
MaxNumSamples = max(sum(indx_Groups));


%Select Pval type
indx_p_Val = strcmpi(p_Val_Type,DATA.ColId);
if any(p_Val_Type)
    p_values = DATA.X(:,indx_p_Val);
else
    error('%s not found',p_Val_Type)
end
switch p_Val_Type
    case 'q t-test'
        p_values =  -log10(p_values);
        p_valueslabel = '-log_1_0(q-value)';
    case 'fdr t-test'
        p_values =  -log10(p_values);
        p_valueslabel = '-log_1_0(fdr p-value)';
    case 'p coxreg DSS'
        p_values =  -log10(p_values);
        p_valueslabel = '-log_1_0(p coxreg DSS)';
    case 'p coxreg PFI'
        p_values =  -log10(p_values);
        p_valueslabel = '-log_1_0(p coxreg PFI)';
end


%Select Y value type
indx_Y_Val = strcmpi(Y_Type,DATA.ColId);
if any(p_Val_Type)
    Y_Val = DATA.X(:,indx_Y_Val);
else
    error('%s not found',Y_Type)
end
switch Y_Type
    case 'Delta Average'
        YLabel = {'\Delta \beta-value'};
    case 'HR coxreg DSS'
        Y_Val =  log2(Y_Val);
        YLabel = {'log_2(HR DSS)'};
    case 'HR coxreg PFI'
        Y_Val =  log2(Y_Val);
        YLabel = {'log_2(HR PFI)'};
end

%
% MaxVal_Y_Significant = max(abs(Y_Val(p_values > p_Val_Cutoff)));
% Y_Val(Y_Val>MaxVal_Y_Significant) = MaxVal_Y_Significant;
% Y_Val(Y_Val<-MaxVal_Y_Significant) = -MaxVal_Y_Significant;


if TopNValues
    p_Val_Cutoff = min(maxk(p_values,TopNValues))
end

if TopPrctile
    p_Val_Cutoff = prctile(p_values,TopPrctile)
end

% get Hyper and Hypo probes
indx_Hyper_S = Y_Val>Y_Val_CutOff & p_values > p_Val_Cutoff;
indx_Hypo_S = Y_Val<-Y_Val_CutOff & p_values > p_Val_Cutoff;

indx_Hyper_pS = Y_Val>0 & p_values > p_Val_Cutoff;
indx_Hypo_pS = Y_Val<0 & p_values > p_Val_Cutoff;



indx_Hyper_ALL = Y_Val>0;
indx_Hypo_ALL = Y_Val<0;
minVal=min(Y_Val(indx_Hypo_ALL));
maxVal=max(Y_Val(indx_Hyper_ALL));
rangeVal = maxVal - minVal;
nudgeVal = rangeVal/5;

nHyper=sum(indx_Hyper_ALL);
nHypo=sum(indx_Hypo_ALL);
nHyper_S = sum(Y_Val>=Y_Val_CutOff & p_values > p_Val_Cutoff);
nHypo_S = sum(Y_Val<=-Y_Val_CutOff & p_values > p_Val_Cutoff);
nHyper_NS = nHyper - nHyper_S;
nHypo_NS = nHypo - nHypo_S;

% Create Plot
cMAP = GetPalette('Tab10');
cMAP = repmat(cMAP,ceil(nGroups/size(cMAP,1)),1);
cMAP = cMAP(1:nGroups,:);

fh = figure('Name','Probe Type Plot','Color','w','Tag','Probe Type Plot','Units','centimeters','GraphicsSmoothing','off');
fWidth = 4+nGroups*1.5;
fHight = 12;
fh.Position(3:4) = [fWidth fHight];
ah = axes(fh,'NextPlot','add','tag','Probe Type Plot','Box','on','FontSize',FontSize,'Linewidth',1,'Units','centimeters');

ah.YLim =[minVal-nudgeVal maxVal+nudgeVal];
ah.OuterPosition(3:4) = [fWidth-2 fHight];
ah.TickLength=[ 0.05/nGroups    0.0];

ShowData=false;
MarkerSize = 25;
Width = 0.45;
JitterWidth = 0.20;
BoxColor = [0.1 0.1 0.1];
TxtStrHyper = cell(3,nGroups);
TxtStrHypo  = cell(3,nGroups);
pHypo       = zeros(nGroups,1);
pHyper      = zeros(nGroups,1);

for i = 1:nGroups
    x_Hypo_ALL = Y_Val(indx_Hypo_ALL & indx_Groups(:,i));
    x_Hypo_S = Y_Val(indx_Hypo_S & indx_Groups(:,i));
    x_Hypo_pS =Y_Val(indx_Hypo_pS & indx_Groups(:,i));
    if Plot_ALL
        Violin({x_Hypo_ALL},i,'ShowData',ShowData,'ViolinColor',{cMAP(i,:)},'ViolinAlpha',{0.5},'Width',Width,'MarkerSize',MarkerSize,'LineWidth', 1,'BoxColor',BoxColor);
    else
        if length(x_Hypo_pS)>0
            Violin({x_Hypo_pS},i,'ShowData',ShowData,'ViolinColor',{cMAP(i,:)},'ViolinAlpha',{0.5},'Width',Width,'MarkerSize',MarkerSize,'LineWidth', 1,'BoxColor',BoxColor);
        end
    end
    x_hypo = i * ones(size(x_Hypo_S)) + (rand(size(x_Hypo_S))-0.5)*(2*JitterWidth);
    scatter(x_hypo,x_Hypo_S,MarkerSize,cMAP(i,:))
    nS=length(x_Hypo_S);
    nALL = length(x_Hypo_ALL);
    [~,pHypo(i),~] = fishertest([nS, nHypo_S-nS; nALL-nS, nHypo_NS-nALL+nS]);
    [stat_str,~] = pval2stars(pHypo(i),[]);

    TxtStrHypo{1,i} =  sprintf('%u',length(x_hypo));
    TxtStrHypo{2,i} =  sprintf('%.2f%%',length(x_hypo)/length(x_Hypo_ALL)*100);
    if ~strcmp('N.S.',stat_str)
        %TxtStr(3,i) = strcat(TxtStr(i),stat_str(i));
        TxtStrHypo(3,i) = stat_str;
    end

    text(i,ah.YLim(1)+0.02,TxtStrHypo(:,i),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',FontSize)


    x_Hyper_ALL = Y_Val(indx_Hyper_ALL & indx_Groups(:,i));
    x_Hyper_S = Y_Val(indx_Hyper_S & indx_Groups(:,i));
    x_Hyper_pS = Y_Val(indx_Hyper_pS & indx_Groups(:,i));
    if Plot_ALL
        Violin({x_Hyper_ALL},i,'ShowData',ShowData,'ViolinColor',{cMAP(i,:)},'ViolinAlpha',{0.5},'Width',Width,'MarkerSize',MarkerSize,'LineWidth', 1,'BoxColor',BoxColor);
    else
        Violin({x_Hyper_pS},i,'ShowData',ShowData,'ViolinColor',{cMAP(i,:)},'ViolinAlpha',{0.5},'Width',Width,'MarkerSize',MarkerSize,'LineWidth', 1,'BoxColor',BoxColor);

    end
    x_hyper = i * ones(size(x_Hyper_S)) + (rand(size(x_Hyper_S))-0.5)*(2*JitterWidth);
    scatter(x_hyper,x_Hyper_S,MarkerSize,cMAP(i,:))
    nS=length(x_Hyper_S);
    nALL = length(x_Hyper_ALL);
    [~,pHyper(i),~] = fishertest([nS, nHyper_S-nS; nALL-nS, nHyper_NS-nALL+nS]);
    [stat_str,~] = pval2stars(pHyper(i),[]);

    TxtStrHyper{1,i} =  sprintf('%u',length(x_hyper));
    TxtStrHyper{2,i} =  sprintf('%.2f%%',length(x_hyper)/length(x_Hyper_ALL)*100);
    if ~strcmp('N.S.',stat_str)
        %TxtStr(3,i) = strcat(TxtStr(i),stat_str(i));
        TxtStrHyper(3,i) = stat_str;
    end

    text(i,ah.YLim(2)-0.02,TxtStrHyper(:,i),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',FontSize)

end

%ah.YTick=ceil((minVal-nudgeVal)*10)/10:0.1:floor((maxVal+nudgeVal)*10)/10;

Prct_Hyper = nHyper_S / nHyper * 100;
Prct_Hypo = nHypo_S / nHypo * 100;

ylabel(ah,YLabel);

line(ah.XLim,[ 0 0 ],'LineWidth',2,'Color',[0 0 0 ],'LineStyle','-');
ah.XLim=[0.5 nGroups+0.5];
ah.XTick=1:nGroups;
ah.XTickLabel=GroupName;
ah.XTickLabelRotation=-45;
ah.XAxis.TickLabelInterpreter='none';
if Y_Val_CutOff
    line(ah.XLim,[ Y_Val_CutOff Y_Val_CutOff],'LineWidth',1,'Color',[1 0 0 ],'LineStyle','--');
    line(ah.XLim,[ -Y_Val_CutOff -Y_Val_CutOff],'LineWidth',1,'Color',[1 0 0 ],'LineStyle','--');
    text(ah.XLim(2)+0.1,Y_Val_CutOff,sprintf('Hyper (%.2f%%)',Prct_Hyper),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',FontSize)
    text(ah.XLim(2)+0.1,-Y_Val_CutOff,sprintf('Hypo (%.2f%%)',Prct_Hypo),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',FontSize)

end


Str = cell(0,0);
[~,cutoff_str] = pval2stars(0.00001,[]);
Str = [{'Fisher’s exact test'}; cutoff_str];
text(ah.XLim(2)+0.1,ah.YLim(2),Str,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',FontSize)

text(ah.XLim(2)+0.1,0,{sprintf('%s < %g',p_valueslabel,p_Val_Cutoff)},'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',FontSize)




