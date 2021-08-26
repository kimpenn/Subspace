function variables = collectOutput( variables, emptyval)
    % Copyright 2016 MathWork, Inc.
    [types,lengths] = getTypesLengths(variables);
    
    block_ends = [0,find(~strcmp(types(1:end-1),types(2:end))),numel(types)];
   
    if numel(types) - nnz(diff(block_ends)>1) < numel(types)
        variables = mergeColumns(variables, types, lengths, block_ends, emptyval);
    end
end

function  variables = mergeColumns(columns, types, lengths, block_ends, emptyval)
    mergedTypes = ~(strcmp(types,{'datetime'})|strcmp(types,{'categorical'}));
    variables = {};
    
    for i = 2:numel(block_ends)
        idx = block_ends(i-1)+1:block_ends(i);

        if mergedTypes(idx(1))
            lengthMismatch = find(abs(diff(lengths(idx)))>0,1);
            if lengthMismatch
                for j = idx(lengthMismatch+1):idx(end)
                    columns{j} = appendVar(columns{j},types{j},emptyval);
                end
            end
            % Concatinate Columns
            variables{end+1} = horzcat(columns{idx}); %#ok<AGROW>
        else
            % datetime and categorical arrays are output as individual columns
            variables(end+1:end+numel(idx)) = columns(idx);
        end
    end
    
end

function [types,lengths] = getTypesLengths(columns)
    types = cell(size(columns));
    lengths = zeros(size(columns));
    
    for i = 1:numel(types)
        var = columns{i};
        types{i} = class(var);
        if strcmp(types{i},'char')
            lengths(i) = size(var,1);
        else
            lengths(i) = numel(var);
        end
    end
end

function var = appendVar(var,type,emptyval)
    switch(type)
        case 'char'
            var(end+1,:) = char(0);
        case 'cell'
            var(end+1,:) = {''};
        case 'string'
            var(end+1,:) = string({''});
        case 'categorical'
            var(end+1,:) = categorical({''});
        case 'datetime'
            var(end+1,:) = NaT;
        otherwise
            var(end+1,:) = emptyval;
    end
end

