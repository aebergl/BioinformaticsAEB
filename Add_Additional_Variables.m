function DATA = Add_Additional_Variables(DATA,DATA2)
% DATA = Add_Additional_Variables(DATA,DATA2)
%
%   Add additional variables to DATA from DATA2



% Check that they are all column vectors
if isrow(DATA.RowId)
    DATA.RowId = DATA.RowId';
end
if isrow(DATA2.RowId)
    DATA2.RowId = DATA2.RowId';
end

if isrow(DATA.ColId)
    DATA.ColId = DATA.ColId';
end
if isrow(DATA2.ColId)
    DATA2.ColId = DATA2.ColId';
end

%Check that the the RowIds are unique
NewCoLId = [DATA.ColId; DATA2.ColId];
if length(NewCoLId) ~= length(unique(NewCoLId))
    warning('DATA and DATA2 have overlapping ColIds. _2 have been added')
    [~,ib] = intersect(DATA.ColId,DATA2.ColId,'stable');
    tmp  = DATA2.ColId;
    tmp(ib) = strcat(tmp(ib),"_2");
    NewCoLId  = [DATA.ColId; tmp];
end

[~,ia,ib] = intersect(DATA.RowId,DATA2.RowId,'stable');

X = ones(DATA.nRow,DATA2.nCol) * NaN;

X(ia,:) = DATA2.X(ib,:);


% Start to Merge the datsets
DATA.X = [DATA.X X];
DATA.nCol = DATA.nCol + DATA2.nCol;
DATA.ColId = NewCoLId;
