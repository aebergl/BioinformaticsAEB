function DATA = CombineVariablesDATA(DATA,VarIds,CombineType,NewVarName)

VarIds = ConvertStr(VarIds,'string');
NewVarName = ConvertStr(NewVarName,'string');
if ~iscolumn(NewVarName)
    NewVarName = NewVarName';
end


DATA_tmp = EditVariablesDATA(DATA,VarIds,'Keep','stable');

switch lower(CombineType)
    case {'+','add','sum'}
        x_new = sum(DATA_tmp.X,2,"omitmissing");
    case {'mean','average'}
        x_new = mean(DATA_tmp.X,2,"omitmissing");
    case 'median'
        x_new = mean(DATA_tmp.X,2,"omitmissing");
    case {'*','prod','multiply'}
        x_new = prod(DATA_tmp.X,2,"omitmissing");
    case {'i','diff'}
        x_new = diff(DATA_tmp.X,2,"omitmissing");
end
nVarAdded = size(x_new,2);
if nVarAdded ~= numel(NewVarName)
    error('There is a mssimatch between number of new variables and New Variables Names given!')

end

DATA.nCol = DATA.nCol + nVarAdded;
DATA.X = [DATA.X x_new];
DATA.ColId = cat(1,DATA.ColId,NewVarName);
if ~isempty(DATA.ColAnnotation)
    DATA.ColAnnotation = append(DATA.ColAnnotation,NewVarName);
end
