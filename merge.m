function [currentPoints, prevAddedIndexes, groups] = merge(index, combined, assignments, maxIndexes, majority, currentPoints, prevAddedIndexes, groups)


for i=1:100
    if any(prevAddedIndexes(:) == maxIndexes(index,i))
        
    else
        currentPoints = vertcat(currentPoints, combined(assignments{maxIndexes(index,i)},:));
        prevAddedIndexes = vertcat(prevAddedIndexes, maxIndexes(index,i));
        groups = vertcat(groups, majority(maxIndexes(index,i)));
        [currentPoints, prevAddedIndexes, groups] = merge(maxIndexes(index,i), combined, assignments, maxIndexes, majority, currentPoints, prevAddedIndexes, groups); 
    end  
end

