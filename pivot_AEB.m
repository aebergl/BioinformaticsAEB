function table_out = pivot_AEB(Group1,Group2)
if iscell(Group1)
        IdInput(cellfun('isempty',Group1)) = {'Empty Cell'};
end
if iscell(Group2)
        IdInput(cellfun('isempty',Group2)) = {'Empty Cell'};
end
% Group2 = strrep(Group2,' - ','_');
% Group2 = strrep(Group2,' ','_');
% Group2 = strrep(Group2,'+','pos');
% Group2 = strrep(Group2,'-','neg');


row_Ids = unique(Group1);
col_Ids = unique(Group2);

X = zeros(length(row_Ids),length(col_Ids));

for i = 1:length(row_Ids)
    indx = strcmp(row_Ids{i},Group1);
    for j = 1:length(col_Ids)
        X(i,j) = sum(strcmp(col_Ids{j},Group2(indx)));
    end
  
end

table_out = array2table(X,'RowNames',row_Ids,'VariableNames',col_Ids);

end

