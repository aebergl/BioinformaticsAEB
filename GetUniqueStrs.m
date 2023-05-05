function [UniqueStr] = GetUniqueStrs(CellData,SplitBy)

% Split all  based on SplitBy
UniqueStr = arrayfun(@(x) split(x,SplitBy,2), CellData, 'UniformOutput', 0);

% Expand to single vector
UniqueStr = [UniqueStr{:}]';

UniqueStr = unique(UniqueStr,'sorted');
% Remove empty cell

if isempty(UniqueStr{1})
    UniqueStr(1) = [];
end
