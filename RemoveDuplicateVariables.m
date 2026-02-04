function DATA = RemoveDuplicateVariables(DATA,VariableId)

VariableId = ConvertStr(VariableId,'char');
VariableIdColumn = strcmpi(VariableId,DATA.ColAnnotationFields);
if ~any(VariableIdColumn)
    error('Given value for Variable Id not found, %s not found in ColAnnotationFields',VariableId)
    return
end


Pattern = ["NA" "---" ""];

DATA =  EditVariablesDATA(DATA,DATA.ColId(matches(DATA.ColAnnotation(:,VariableIdColumn),Pattern)),'remove');


[IDs,numID]  = GroupCount(DATA.ColAnnotation(:,VariableIdColumn),0);

MultipleIDs = find(numID>1);
for i=1:length(MultipleIDs)
    id = IDs(MultipleIDs(i));
    DATA_tmp = EditVariablesDATA(DATA,id,'Keep','VariableIdentifier',VariableId);
    [~,indx] = max(mean(DATA_tmp.X,1,'omitnan'));
    ColIdsToRemove = DATA_tmp.ColId(~matches(DATA_tmp.ColId,DATA_tmp.ColId(indx)));
    DATA = EditVariablesDATA(DATA,ColIdsToRemove,'Remove');
end

