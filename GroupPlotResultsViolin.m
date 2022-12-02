function fh = GroupPlotResultsViolin(DATA,GroupId,Groups,p_Val_Cutoff,p_Val_Type,DeltaBeta_Cutoff)
FontSize = 10;
Plot_ALL=false;

if ~iscell(Groups)
    Groups = {Groups};
end



nVals = DATA.nRow;

GroupColumn = strcmp(GroupId,DATA.RowAnnotationFields);
if ~any(GroupColumn)
    error('No group Column Found')
end
nGroups = numel(Groups);



numValuesPerGroup = zeros(nGroups,1);
indx_Groups=false(DATA.nRow,nGroups);
GroupVariable = cell(DATA.nRow,1);
GroupName = cell(numGroups1,1);

for i=1:nGroups
    tmp_name = Groups{i};
    if ~iscell(tmp_name)
        tmp_name = {tmp_name};
    end
    SampleIndxMatrix = false(DATA.NumSamples,numel(tmp_name));
    for j=1:length(tmp_name)
        SampleIndxMatrix(:,j) = strcmp(tmp_name{j},DATA.RowAnnotation(:,GroupColumn));
    end
    indx_Groups(:,i) = any(SampleIndxMatrix,2);
    GroupVariable(indx_Groups(:,i)) = tmp_name(1);
    numValuesPerGroup(i) = sum(indx_Groups(:,i));
    GroupName(i) = tmp_name(1);
end
MaxNumSamples = max(sum(indx_Groups));


indx_Groups = false(nVals,nGroups);




if strcmpi('q',p_Val_Type)
    p_values = DATA.RES(nRES).Q_TT;
else
    p_values = DATA.RES(nRES).pTT;
end



% get Hyper and Hypo probes
indx_Hyper_S = DATA.RES(nRES).MeanDiff>DeltaBeta_Cutoff & p_values < p_Val_Cutoff;
indx_Hypo_S = DATA.RES(nRES).MeanDiff<-DeltaBeta_Cutoff & p_values < p_Val_Cutoff;

indx_Hyper_pS = DATA.RES(nRES).MeanDiff>0 & p_values < p_Val_Cutoff;
indx_Hypo_pS = DATA.RES(nRES).MeanDiff<0 & p_values < p_Val_Cutoff;



indx_Hyper_ALL = DATA.RES(nRES).MeanDiff>0;
indx_Hypo_ALL = DATA.RES(nRES).MeanDiff<0;
minVal=min(DATA.RES(nRES).MeanDiff(indx_Hypo_ALL));
maxVal=max(DATA.RES(nRES).MeanDiff(indx_Hyper_ALL));
rangeVal = maxVal - minVal;
nudgeVal = rangeVal/5;

nHyper=sum(indx_Hyper_ALL);
nHypo=sum(indx_Hypo_ALL);
nHyper_S = sum(DATA.RES(nRES).MeanDiff>=DeltaBeta_Cutoff & p_values < p_Val_Cutoff);
nHypo_S = sum(DATA.RES(nRES).MeanDiff<=-DeltaBeta_Cutoff & p_values < p_Val_Cutoff);
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
JitterWidth = 0.10;
BoxColor = [0.1 0.1 0.1];
TxtStrHyper = cell(3,nGroups);
TxtStrHypo  = cell(3,nGroups);
pHypo       = zeros(nGroups,1);
pHyper      = zeros(nGroups,1);
for i = 1:nGroups
    x_Hypo_ALL = DATA.RES(nRES).MeanDiff(indx_Hypo_ALL & indx_Groups(:,i));
    x_Hypo_S = DATA.RES(nRES).MeanDiff(indx_Hypo_S & indx_Groups(:,i));
    x_Hypo_pS = DATA.RES(nRES).MeanDiff(indx_Hypo_pS & indx_Groups(:,i));
    if Plot_ALL
        Violin({x_Hypo_ALL},i,'ShowData',ShowData,'ViolinColor',{cMAP(i,:)},'ViolinAlpha',{0.5},'Width',Width,'MarkerSize',MarkerSize,'LineWidth', 1,'BoxColor',BoxColor);
    else
        Violin({x_Hypo_pS},i,'ShowData',ShowData,'ViolinColor',{cMAP(i,:)},'ViolinAlpha',{0.5},'Width',Width,'MarkerSize',MarkerSize,'LineWidth', 1,'BoxColor',BoxColor);
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
    
    
    x_Hyper_ALL = DATA.RES(nRES).MeanDiff(indx_Hyper_ALL & indx_Groups(:,i));
    x_Hyper_S = DATA.RES(nRES).MeanDiff(indx_Hyper_S & indx_Groups(:,i));
    x_Hyper_pS = DATA.RES(nRES).MeanDiff(indx_Hyper_pS & indx_Groups(:,i));
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

ah.YTick=ceil((minVal-nudgeVal)*10)/10:0.1:floor((maxVal+nudgeVal)*10)/10;



ylabel(ah,'\Delta \beta-value');

line(ah.XLim,[ DeltaBeta_Cutoff DeltaBeta_Cutoff],'LineWidth',1,'Color',[1 0 0 ],'LineStyle','--');
line(ah.XLim,[ -DeltaBeta_Cutoff -DeltaBeta_Cutoff],'LineWidth',1,'Color',[1 0 0 ],'LineStyle','--');

line(ah.XLim,[ 0 0 ],'LineWidth',2,'Color',[0 0 0 ],'LineStyle','-');
ah.XLim=[0.5 nGroups+0.5];
ah.XTick=1:nGroups;
ah.XTickLabel=Groups;
ah.XTickLabelRotation=-45;
ah.XAxis.TickLabelInterpreter='none';

        Prct_Hyper = nHyper_S / nHyper * 100;
        Prct_Hypo = nHypo_S / nHypo * 100;

        Str = cell(0,0);
            [~,cutoff_str] = pval2stars(0.00001,[]);
        Str = [{'Fisher’s exact test'}; cutoff_str];
        text(ah.XLim(2)+0.1,ah.YLim(2),Str,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',FontSize)

        text(ah.XLim(2)+0.1,0,{sprintf('%s < %g',p_Val_Type,p_Val_Cutoff)},'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',FontSize)

        text(ah.XLim(2)+0.1,DeltaBeta_Cutoff,sprintf('Hyper (%.2f%%)',Prct_Hyper),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',FontSize)
        text(ah.XLim(2)+0.1,-DeltaBeta_Cutoff,sprintf('Hypo (%.2f%%)',Prct_Hypo),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',FontSize)



