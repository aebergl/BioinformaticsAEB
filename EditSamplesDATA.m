function DATA  = EditSamplesDATA(DATA,InputIds,KeepRemove)

IdUnique = unique(InputIds);
if length(IdUnique) < length(InputIds)
    warning('WARNING! Not all input IDs are unique in EditSamplesDATA')
end

indx = ismember(DATA.RowId,InputIds);

if any(indx)
    switch lower(KeepRemove)
        case 'keep'
            DATA.X = DATA.X(indx,:);
            DATA.RowId = DATA.RowId(indx);
            DATA.nRow = size(DATA.X,1);
            DATA.RowAnnotation = DATA.RowAnnotation(indx,:);
            if isfield(DATA,'SURVIVAL')
                DATA.SURVIVAL.RowId =  DATA.SURVIVAL.RowId(indx);
                DATA.SURVIVAL.SurvEvent =  DATA.SURVIVAL.SurvEvent(indx,:);
                DATA.SURVIVAL.SurvTime =  DATA.SURVIVAL.SurvTime(indx,:);
            end
        case 'remove'
            %Remove samples
            DATA.X(indx,:) = [];
            DATA.RowId(indx) = [];
            DATA.nRow = size(DATA.X,1);
            DATA.RowAnnotation(indx,:) = [];
            if isfield(DATA,'SURVIVAL')
                DATA.SURVIVAL.RowId(indx) = [];
                DATA.SURVIVAL.SurvEvent(indx,:) =  [];
                DATA.SURVIVAL.SurvTime(indx,:) =  [];
            end
    end
end

end