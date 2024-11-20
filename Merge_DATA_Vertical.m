function DATA = Merge_DATA_Vertical(DATA,DATA2,DataSetId,DataSetNames)
% DATA = Merge_DATA_Vertical(DATA,DATA2)
%
%   Merge two data strcucures vertically


% Check for correct size
if DATA.nCol ~= DATA2.nCol
    error('The two datasets have different number of columns')
end

% Check for matching ColId
if ~all(string(DATA.ColId) == string(DATA2.ColId))
    error('DATA and DATA2 have different ColIds')
end

% Check for matching RowAnnotationFields
if ~all(string(DATA.RowAnnotationFields) == string(DATA2.RowAnnotationFields))
    error('DATA and DATA2 have different RowAnnotationFields')
end

% Check that they are all column vectors
if isrow(DATA.RowId)
    DATA.RowId = DATA.RowId';
end
if isrow(DATA2.RowId)
    DATA2.RowId = DATA2.RowId';
end

%Check that the the RowIds are unique
NewRowId = [DATA.RowId; DATA2.RowId];
if length(NewRowId) ~= length(unique(NewRowId))
    error('DATA and DATA2 have overlapping RowIds. They need to be unique')
end

% Start to Merge the datsets
DATA.X = [DATA.X; DATA2.X];
DATA.nRow = DATA.nRow + DATA2.nRow;
DATA.RowId = NewRowId;
DATA.RowAnnotation = [DATA.RowAnnotation; DATA2.RowAnnotation];