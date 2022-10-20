function [IDs,numID] = GroupCount(IdInput,PrintResults)

if nargin < 2
    PrintResults = 0;
end
if iscell(IdInput)
        IdInput(cellfun('isempty',IdInput)) = {'Empty Cell'};
else

end
    

IDs = unique(IdInput,'stable');
numID = zeros(length(IDs),1);

for i = 1:length(IDs)
    if iscell(IdInput)
        numID(i) = sum(strcmp(IDs{i},IdInput));
    else
        numID(i) = sum(IdInput == IDs(i));
    end

    if PrintResults
         if iscell(IdInput)
        fprintf('%s:\t%u\n',IDs{i},numID(i))
         else
              fprintf('%i:\t%u\n',IDs(i),numID(i))
         end
    end
end
