function DATA_Out = PCAmodel2DATA(PCAmodel,TypeFlag,DATA)

switch lower(TypeFlag)
    case 'scores'
        [nRow, nCol]= size(PCAmodel.T);
        DATA_Out = CreateDataStructure(nRow,nCol,1,1);
        DATA_Out.X = PCAmodel.T;
        for i=1:length(PCAmodel.ExplVar)
            DATA_Out.ColId(i) = sprintf('PC%u (%.1f%%)',i,PCAmodel.ExplVar(i));
        end
        DATA_Out.RowId = string(1:nRow)';
        if ~isempty(DATA) % Merge in DATA
            if isempty(PCAmodel.row_rem)
                PCAmodel.row_rem = zeros(nRow,1);
            end
            if length(PCAmodel.row_rem) == DATA.nRow
                DATA_Out.RowId = DATA.RowId(~PCAmodel.row_rem);
                if ~isempty(DATA.RowAnnotationFields)
                    DATA_Out.RowAnnotationFields = DATA.RowAnnotationFields;
                    DATA_Out.RowAnnotation = DATA.RowAnnotation(~PCAmodel.row_rem,:);
                end
                if isfield(DATA,'SURVIVAL')
                    DATA_Out.SURVIVAL.RowId =  DATA.SURVIVAL.RowId(~PCAmodel.row_rem);
                    DATA_Out.SURVIVAL.SurvEvent =  DATA.SURVIVAL.SurvEvent(~PCAmodel.row_rem,:);
                    DATA_Out.SURVIVAL.SurvTime =  DATA.SURVIVAL.SurvTime(~PCAmodel.row_rem,:);
                end
            else
                error("Not maching number of rows")
            end
        end
    case 'loadings'
        [nRow, nCol]= size(PCAmodel.P);
        DATA_Out = CreateDataStructure(nRow,nCol,1,1);
        DATA_Out.X = PCAmodel.P;
        for i=1:length(PCAmodel.ExplVar)
            DATA_Out.ColId(i) = sprintf('PC Loading %u (%.1f%%)',i,PCAmodel.ExplVar(i));
        end
        DATA_Out.RowId = string(1:nRow)';
        if ~isempty(DATA) % Merge in DATA
            if isempty(PCAmodel.row_rem)
                PCAmodel.col_rem = zeros(nRow,1);
            end

            if length(PCAmodel.col_rem) == DATA.nCol
                DATA_Out.RowId = DATA.ColId(~PCAmodel.col_rem);
                if ~isempty(DATA.ColAnnotationFields)
                    DATA_Out.RowAnnotationFields = DATA.ColAnnotationFields;
                    DATA_Out.RowAnnotation = DATA.ColAnnotation(~PCAmodel.col_rem,:);
                end
            else
                error("Not maching number of rows")
            end
        end
    otherwise
        error('TypeFlag must be "Scores" or "Loadings"')
end