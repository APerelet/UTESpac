function fullStruct = structConcat(struct1, struct2)
% vertically concatnates all fields within 2 structures.  Struct 2 will be concatnated to the end of structg 1.Field
% names and field columns must be perfectly consistent accross structures.  Cells are not vertically concatnated!

% struct1 = ES5F;
% struct2 = ES5S;

% initialize concatnated structure
fullStruct = struct;

% get all fields
allFields = fields(struct1);

for i = 1:numel(allFields)
    
    % check for cell
    if iscell(struct1.(allFields{i}))
        % if cell, take field from struct1
        fullStruct.(allFields{i}) = struct1.(allFields{i});
        continue
    end
    
    % vertically concatnate fields
    fullStruct.(allFields{i}) = [struct1.(allFields{i}); struct2.(allFields{i})];
end