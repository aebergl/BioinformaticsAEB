function fh = Plot_GSEA_Heatmap(DATA,XId,YId,ColorId,PvalId)


Pval_cutoff = 0.01;

indx_Pval = contains(DATA.ColId,PvalId);

Pval = DATA.X(:,indx_Pval);

indx_Yid = contains(DATA.RowAnnotationFields,YId);

SignificantGeneSets = unique((DATA.RowAnnotation(Pval<Pval_cutoff,indx_Yid)));

DATA = EditSamplesDATA(DATA,SignificantGeneSets,'Keep','SampleIdentifier',YId,'Stable');
Pval = DATA.X(:,indx_Pval);
indx_ColorId = contains(DATA.ColId,ColorId);
%DATA.X(Pval>Pval_cutoff,indx_ColorId) = NaN;


% IDs =  ["metabolic" "proliferation" "signaling" "cellular component" "development" "DNA damage" "immune" "pathway"];
% DATA = EditSamplesDATA(DATA,IDs,'Keep','SampleIdentifier','Group','Stable');

ListOfPathways=unique(DATA.RowAnnotation(:,indx_Yid),'stable');


T1=array2table(DATA.RowAnnotation,RowNames=DATA.RowId,VariableNames=DATA.RowAnnotationFields);
T2=array2table(DATA.X,RowNames=DATA.RowId,VariableNames=DATA.ColId);
T=join(T1,T2,'Keys', 'Row');

fh = figure('Name','Probe correlation heatmap','Color','w','Tag','Probe correlation heatmap');



 h=heatmap(fh,T,XId,YId,'ColorVariable',ColorId,'Colormap',colorcet('D1'),'ColorLimits',[-3 3],...
     'MissingDataColor','w','YDisplayData',ListOfPathways);

h.MissingDataLabel = 'N.S.';
h.Title = [];
h.XLabel=[];
h.YLabel=[];