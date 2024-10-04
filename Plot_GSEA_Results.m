function fh = Plot_GSEA_Results(File_pos,File_neg,X_value_name,Size_value_Name,CollectionType,PvalCutOff,Nmax,fWidth,fHight)

FontSize = 7;
minSize = 20;
maxSize = 100;
LineWidth = 0.5;
GridLines = 'on';
RightMargin = 0.8;

NameField = "NAME";
SizeCutOffs =[1 0.05 0.01 0.001 0.0001];
LegendSizeVal = [5 20 40 60 80];
% SizeCutOffs =[1 0.2 0.1 0.05 ];
% LegendSizeVal = [5 40 60 100];
% SizeCutOffs =[1 0.15 0.1 0.05 0.01 ];
% LegendSizeVal = [5 20 40 60 80];

FlipNeg = true;
ColorBy = "Enrichment";
ColorBy = "Type";

xVal_pos = [];
SizeVal_pos = [];
YtickLabelTxt_pos = [];

xVal_neg = [];
SizeVal_neg = [];
YtickLabelTxt_neg = [];

% Read Files
%Read input data
if ~ isempty(File_pos)
    opts_pos = detectImportOptions(File_pos,'FileType','delimitedtext','VariableNamingRule','preserve',"TextType","string","Delimiter","\t");
    T_pos = readtable(File_pos,opts_pos);
    xVal_pos = table2array(T_pos(:,X_value_name));
    SizeVal_pos = table2array(T_pos(:,Size_value_Name));
    YtickLabelTxt_pos = table2array(T_pos(:,NameField));
    indx = SizeVal_pos < PvalCutOff;
    xVal_pos = xVal_pos(indx);
    SizeVal_pos = SizeVal_pos(indx);
    YtickLabelTxt_pos = YtickLabelTxt_pos(indx);
    if length(xVal_pos) > Nmax
        xVal_pos = xVal_pos(1:Nmax);
        SizeVal_pos = SizeVal_pos(1:Nmax);
        YtickLabelTxt_pos = YtickLabelTxt_pos(1:Nmax);
    end

end
nPos = length(xVal_pos);

if ~ isempty(File_neg)
    opts_neg=detectImportOptions(File_neg,'FileType','delimitedtext','VariableNamingRule','preserve',"TextType","string","Delimiter","\t");
    T_neg=readtable(File_neg,opts_neg);
    xVal_neg = table2array(T_neg(:,X_value_name));
    SizeVal_neg = table2array(T_neg(:,Size_value_Name));
    YtickLabelTxt_neg = table2array(T_neg(:,NameField));
    indx = SizeVal_neg < PvalCutOff;
    xVal_neg = xVal_neg(indx);
    SizeVal_neg = SizeVal_neg(indx);
    YtickLabelTxt_neg = YtickLabelTxt_neg(indx);

    if length(xVal_neg) > Nmax
        xVal_neg = xVal_neg(1:Nmax);
        SizeVal_neg = SizeVal_neg(1:Nmax);
        YtickLabelTxt_neg = YtickLabelTxt_neg(1:Nmax);
    end
    if FlipNeg
        xVal_neg = flipud(xVal_neg);
        SizeVal_neg = flipud(SizeVal_neg);
        YtickLabelTxt_neg = flipud(YtickLabelTxt_neg);

    end
end
nNeg = length(xVal_neg);

% Combine pos & neg
YtickLabelTxt = [YtickLabelTxt_pos; YtickLabelTxt_neg];
xVal = [xVal_pos; xVal_neg];
SizeVal = [SizeVal_pos; SizeVal_neg];
nGroups = nPos + nNeg;

% Fix Gene set names
switch lower(CollectionType)
    case 'hallmark'
        GSEA_Name ={
            'HALLMARK_ADIPOGENESIS','Adipogenesis','development';
            'HALLMARK_ALLOGRAFT_REJECTION','Allograft Rejection','immune';
            'HALLMARK_ANDROGEN_RESPONSE','Androgen Response','signaling';
            'HALLMARK_ANGIOGENESIS','Angiogenesis','development';
            'HALLMARK_APICAL_JUNCTION','Apical Junction','cellular component';
            'HALLMARK_APICAL_SURFACE','Apical Surface','cellular component';
            'HALLMARK_APOPTOSIS','Apoptosis','pathway';
            'HALLMARK_BILE_ACID_METABOLISM','Bile Acid Metabolism','metabolic';
            'HALLMARK_CHOLESTEROL_HOMEOSTASIS','Cholesterol Homeostasis','metabolic';
            'HALLMARK_COAGULATION','Coagulation','immune';
            'HALLMARK_COMPLEMENT','Complement','immune';
            'HALLMARK_DNA_REPAIR','DNA Repair','DNA damage';
            'HALLMARK_E2F_TARGETS','E2F Targets','proliferation';
            'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','EMT','development';
            'HALLMARK_ESTROGEN_RESPONSE_EARLY','Estrogen Response Early','signaling';
            'HALLMARK_ESTROGEN_RESPONSE_LATE','Estrogen Response Late','signaling';
            'HALLMARK_FATTY_ACID_METABOLISM','Fatty Acid Metabolism','metabolic';
            'HALLMARK_G2M_CHECKPOINT','G2M Checkpoint','proliferation';
            'HALLMARK_GLYCOLYSIS','Glycolysis','metabolic';
            'HALLMARK_HEDGEHOG_SIGNALING','Hedgehog Signaling','signaling';
            'HALLMARK_HEME_METABOLISM','Heme Metabolism','metabolic';
            'HALLMARK_HYPOXIA','Hypoxia','pathway';
            'HALLMARK_IL2_STAT5_SIGNALING','IL2 STAT5 Signaling','signaling';
            'HALLMARK_IL6_JAK_STAT3_SIGNALING','IL6 JAK STAT3 Signaling','immune';
            'HALLMARK_INFLAMMATORY_RESPONSE','Inflammatory Response','immune';
            'HALLMARK_INTERFERON_ALPHA_RESPONSE','Interferon alpha Response','immune';
            'HALLMARK_INTERFERON_GAMMA_RESPONSE','Interferon gamma Response','immune';
            'HALLMARK_KRAS_SIGNALING_DN','KRAS Signaling Dn','signaling';
            'HALLMARK_KRAS_SIGNALING_UP','KRAS Signaling Up','signaling';
            'HALLMARK_MITOTIC_SPINDLE','Mitotic Spindle','proliferation';
            'HALLMARK_MTORC1_SIGNALING','MTORC1 Signaling','signaling';
            'HALLMARK_MYC_TARGETS_V1','MYC Targets V1','proliferation';
            'HALLMARK_MYC_TARGETS_V2','MYC Targets V2','proliferation';
            'HALLMARK_MYOGENESIS','Myogenesis','development';
            'HALLMARK_NOTCH_SIGNALING','NOTCH Signaling','signaling';
            'HALLMARK_OXIDATIVE_PHOSPHORYLATION','Oxidative Phosphorylation','metabolic';
            'HALLMARK_P53_PATHWAY','P53 Pathway','proliferation';
            'HALLMARK_PANCREAS_BETA_CELLS','Pancreas beta Cell','development';
            'HALLMARK_PEROXISOME','Peroxisome','cellular component';
            'HALLMARK_PI3K_AKT_MTOR_SIGNALING','PI3K AKT MTOR Signaling','signaling';
            'HALLMARK_PROTEIN_SECRETION','Protein Secretion','pathway';
            'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY','Reactive Oxygen Species','pathway';
            'HALLMARK_SPERMATOGENESIS','Spermatogenesis','development';
            'HALLMARK_TGF_BETA_SIGNALING','TGF beta Signaling','signaling';
            'HALLMARK_TNFA_SIGNALING_VIA_NFKB','TNFa Signaling Via NFkB','signaling';
            'HALLMARK_UNFOLDED_PROTEIN_RESPONSE','Unfolded Protein Response','pathway';
            'HALLMARK_UV_RESPONSE_DN','UV Response Down','DNA damage';
            'HALLMARK_UV_RESPONSE_UP','UV Response Up','DNA damage';
            'HALLMARK_WNT_BETA_CATENIN_SIGNALING','WNT beta Catenin Signaling','signaling';
            'HALLMARK_XENOBIOTIC_METABOLISM','Xenobiotic Metabolism','metabolic';
            };
        [~,indx_new]=ismember(YtickLabelTxt,GSEA_Name(:,1));

        YtickLabelTxt = GSEA_Name(indx_new,2);
        ColorGroups = GSEA_Name(indx_new,3);
        AllTypes = {'cellular component','development','DNA damage','immune','metabolic','pathway','proliferation','signaling'};
        % YtickLabelTxt = strrep(YtickLabelTxt,'HALLMARK_','');
        % YtickLabelTxt = strrep(YtickLabelTxt,'_',' ');
        % YtickLabelTxt = lower(YtickLabelTxt);
        % YtickLabelTxt = upper(extractBefore(YtickLabelTxt,2)) + extractAfter(YtickLabelTxt,1);
    case 'reactome'
        YtickLabelTxt = strrep(YtickLabelTxt,'REACTOME_','');
        YtickLabelTxt = strrep(YtickLabelTxt,'_',' ');
        YtickLabelTxt = lower(YtickLabelTxt);
        YtickLabelTxt = upper(extractBefore(YtickLabelTxt,2)) + extractAfter(YtickLabelTxt,1);
end


% SizeValPlot = rescale(SizeVal,minSize,maxSize);
% LegendSizeValPlot = rescale(LegendSizeVal,minSize,maxSize);
% Color type

switch lower(ColorBy)

    case 'enrichment'
        ColorVal = [1*ones(nPos,1); 2*ones(nNeg,1)];
        Cmap = GetPalette('Lancet', [2 1]);
        ColorVal = Cmap(ColorVal,:);
    case 'type'
        ColorVal = categorical(ColorGroups,AllTypes);
        ColorVal = double(ColorVal);
        %Cmap = GetPalette('aeb01',[ 7 6 3 4 8 9 11 5]);
        Cmap = GetPalette('Tab10',[ 7 3 5 10 9 8 2 6]);
        ColorVal = Cmap(ColorVal,:);
        YtickLabelTxt(1:nPos) = strcat('\color{red}',YtickLabelTxt(1:nPos));
        YtickLabelTxt(nPos+1:end) = strcat('\color{blue}',YtickLabelTxt(nPos+1:end));
end
% Assign size for each sample
LegendSizeValMat = sum(repmat(SizeCutOffs,nGroups,1) >= repmat(SizeVal,1,length(SizeCutOffs)),2);
SizeValPlot = LegendSizeVal(LegendSizeValMat);

minVal=min(xVal);
maxVal=max(xVal);
rangeVal = maxVal - minVal;
nudgeVal = rangeVal/10;


fh=figure('Name','GSEA Plot','Color','w','Tag','GSEA Plot figure','Units','inches');
fh.Position(3:4) = [fWidth fHight];
ah = axes(fh,'NextPlot','add','tag','Volcano Plot','box','on','Layer','top','FontSize',FontSize,'Units','inches',...
    'PositionConstraint','outerposition','Clipping','off','YDir','reverse');

ah.LineWidth = LineWidth;
ah.XGrid = GridLines;
ah.YGrid = GridLines;
ah.YDir = 'Reverse';
ah.YTick = 1:nGroups;
ah.YLim = [0.5 nGroups+0.5];
ah.YTickLabel = YtickLabelTxt;
ah.Colormap = Cmap;
ah.ColorOrder = Cmap;
%ah.CLim = CLim;
ah.XLim =[minVal-nudgeVal maxVal+nudgeVal];
ah.OuterPosition(3:4) = [fWidth-RightMargin fHight];
ah.TickLength=[ 0.05/nGroups    0.0];

switch lower(X_value_name)
    case 'es'
        XlabelTxt = "Enrichment score (ES)";
    case 'nes'
        XlabelTxt = "Normalized enrichment score (NES)";
end
xlabel(ah,XlabelTxt);

if nPos && nNeg
    line(ah,[0 0],ah.YLim,'Linewidth',0.5,'Color','k','LineStyle','-')
end


sh = scatter(xVal,1:nGroups,SizeValPlot,ColorVal,'MarkerFaceColor','flat','MarkerEdgeColor',[0.1 0.1 0.1]);

StepVal = 1.1;

% Add color legend
YPos_Size = 0;
YPos_Text = 0;
YPos_Color = 0;

YPos_Text = 1:StepVal:StepVal*3;


text(ah.XLim(2)+1.8*nudgeVal,YPos_Text(1),CollectionType,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',FontSize,'FontWeight','bold')
text(ah.XLim(2)+1.8*nudgeVal,YPos_Text(2),'Up-regulated','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',FontSize,'Color','red')
text(ah.XLim(2)+1.8*nudgeVal,YPos_Text(3),'Down-regulated','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',FontSize,'Color','blue')


switch lower(ColorBy)

    case 'enrichment'
        text(ah.XLim(2)+1.4*nudgeVal,YPos_Size(end)+1+(2*StepVal),'Enrichment','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',FontSize)
        scatter(ah,ah.XLim(2)+nudgeVal,YPos_Size(end)+1+(3*StepVal),LegendSizeVal(end),Cmap(1,:),'MarkerFaceColor','flat','MarkerEdgeColor',[0.1 0.1 0.1]);
        text(ah.XLim(2)+1.6*nudgeVal,YPos_Size(end)+1+(3*StepVal),'Up-regulated','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',FontSize)
        scatter(ah,ah.XLim(2)+nudgeVal,YPos_Size(end)+1+(4.2*StepVal),LegendSizeVal(end),Cmap(2,:),'MarkerFaceColor','flat','MarkerEdgeColor',[0.1 0.1 0.1]);
        text(ah.XLim(2)+1.6*nudgeVal,YPos_Size(end)+1+(4.2*StepVal),'Down-regulated','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',FontSize)
    case 'type'
        ColorGroupslegend = unique(ColorGroups,'stable');
        ColorValLegend = unique(ColorVal,'stable','rows');
%        ColorValLegend = double(ColorValLegend);
        YPos_Color = 1:StepVal:StepVal*length(ColorValLegend);
        YPos_Color = YPos_Color + YPos_Text(end) + 1;
        for i=1:length(ColorValLegend)
            scatter(ah,ah.XLim(2)+nudgeVal*0.9,YPos_Color(i),LegendSizeVal(end-1),ColorValLegend(i,:),'MarkerFaceColor','flat','MarkerEdgeColor',[0.1 0.1 0.1]);
            text(ah.XLim(2)+1.3*nudgeVal,YPos_Color(i),ColorGroupslegend{i},'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',FontSize)
        end
end

% Add size legend
YPos_Size = 1:StepVal:StepVal*length(LegendSizeVal);
YPos_Size = YPos_Size + YPos_Color(end) + 3;
shl = scatter(ah,ah.XLim(2)+nudgeVal,YPos_Size,LegendSizeVal,[0 0 0]);
text(ah.XLim(2)+1.3*nudgeVal,YPos_Size(1)-StepVal,Size_value_Name,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',FontSize)
for i = 1:length(LegendSizeVal)
    if i==1
        txt_str = 'N.S.';
    elseif i==length(LegendSizeVal)
        txt_str = strcat("<",num2str(SizeCutOffs(i)));
    else
        txt_str = num2str(SizeCutOffs(i));
    end
    text(ah.XLim(2)+1.6*nudgeVal,YPos_Size(i),txt_str,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',FontSize)
end

