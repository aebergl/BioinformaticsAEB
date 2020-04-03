function DATA  = EditSamplesDATA(DATA,InputIds,KeepRemove)

IdUnique = unique(InputIds);
if length(IdUnique) < length(InputIds)
    fprintf('WARNING! Not all input IDs are unique in EditSamplesDATA\n')
end

indx = ismember(DATA.RowId,InputIds);

if any(indx)
    disp(' ')
    disp('***   Error   ***')
    disp('No matching ids found')
    disp(' ')
end

if strcmpi('Keep',KeepRemove)
    DATA.X = DATA.X(DATA_indx,:);
    DATA.SampleId = DATA.SampleId(DATA_indx);
    DATA.NumSamples = size(DATA.X,1);
    DATA.SampleAnnotation = DATA.SampleAnnotation(DATA_indx,:);
    if isfield(DATA,'SURVIVAL')
        DATA.SURVIVAL.time =  DATA.SURVIVAL.time(DATA_indx,:);
        DATA.SURVIVAL.Status =  DATA.SURVIVAL.Status(DATA_indx,:);
    end
elseif strcmpi('Remove',KeepRemove)
    %Remove samples
    DATA.X(DATA_indx,:) = [];
    DATA.SampleId(DATA_indx) = [];
    DATA.NumSamples = size(DATA.X,1);
    DATA.SampleAnnotation(DATA_indx,:) = [];
    if isfield(DATA,'SURVIVAL')
        DATA.SURVIVAL.time(DATA_indx,:) =  [];
        DATA.SURVIVAL.Status(DATA_indx,:) =  [];
    end

else
    disp(' ')
    disp('Second variable needs to be Keep or Remove')
    disp('exiting!!!')
    disp(' ')
    
    return
end


    
    
end