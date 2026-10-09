function [DATA_scores, DATA_stats]  = ScoreSignatureFromDATA(DATA,VariableId,SIGs,SigId)

ScaleXflag = true;
nVarCutOff = 5;

nSig = numel(SIGs);

DATA_scores = CreateDataStructure(DATA.nRow,nSig,[],[]);
DATA_scores.RowId = DATA.RowId;
DATA_scores.RowAnnotationFields = DATA.RowAnnotationFields;
DATA_scores.RowAnnotation = DATA.RowAnnotation;
DATA_scores.ColAnnotationFields = ["Description", "Source"];
DATA_scores.ColAnnotation = strings(nSig,2);


% Create Stats Strucure
StatsVar = ["%PC1" "PC1/PC2" "nVar"];

DATA_stats = CreateDataStructure(nSig,length(StatsVar),[],[]);



nSigsCounter = 0;
for i=1:nSig
    indx = strcmp(SigId,SIGs(i).Identifiers);
    if any(indx)
        ids = SIGs(i).IDs(:,indx);
        indx_empty = (ids == "");
        ids(indx_empty) = [];
        ids = unique(ids);
        if length(ids) >= nVarCutOff
            DATA_tmp = EditVariablesDATA(DATA,ids,'Keep','VariableIdentifier',VariableId,'Stable');
            if ~isempty(DATA_tmp) && DATA_tmp.nCol >= nVarCutOff
                nSigsCounter = nSigsCounter + 1;
                PCA_model = NIPALS_PCA(DATA_tmp.X,'NumComp',2,'ScaleX',ScaleXflag);
                DATA_scores.X(:,nSigsCounter) = PCA_model.T(:,1);
                DATA_scores.ColId(nSigsCounter) = SIGs(i).Name;
                DATA_scores.ColAnnotation(nSigsCounter,1) = SIGs(i).Description;
                DATA_scores.ColAnnotation(nSigsCounter,2) = SIGs(i).Source;

                DATA_stats.X(nSigsCounter,1) = PCA_model.ExplVar(1);
                DATA_stats.X(nSigsCounter,2) = PCA_model.ExplVar(1) / PCA_model.ExplVar(2);
                DATA_stats.X(nSigsCounter,3) = DATA_tmp.nCol;
                DATA_stats.RowId(nSigsCounter) = SIGs(i).Name;
            end
        end
    end
end
DATA_scores.nCol = nSigsCounter;
DATA_scores.X = DATA_scores.X(:,1:nSigsCounter);
DATA_scores.ColId = DATA_scores.ColId(1:nSigsCounter);
DATA_scores.ColAnnotation = DATA_scores.ColAnnotation(1:nSigsCounter,:);

DATA_stats.nRow = nSigsCounter;
DATA_stats.X = DATA_stats.X(1:nSigsCounter,:);
DATA_stats.RowId = DATA_stats.RowId(1:nSigsCounter);
