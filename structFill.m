function filledStrct = structFill(strct,fillDates,avgPer)

% find all fields
allFields = fields(strct);

% initialize output structure
filledStrct = struct;

% find first logger table
firstLoggerTableRow = find(~cellfun(@isempty,strfind(allFields,'eader')),1,'first')-1;

% find number of rows per day
numRows = sum(floor(strct.(allFields{firstLoggerTableRow})(:,1))==floor(strct.(allFields{firstLoggerTableRow})(1,1)))+1;

% iterate through all fields
for i = 1:numel(allFields)
    
    % find local variable
    localVar = strct.(allFields{i});
    
    % if localVar is a cell, store in filledStrct and continue
    if iscell(localVar)
        filledStrct.(allFields{i}) = localVar;
        continue
    end
    
    % find number of columns in localVar
    numCols = size(localVar,2);
    
    % check to make sure timestamp exists in col 1, if not, find appropriate time stamp
    if localVar(1,1) > datenum(2000,0,0) && localVar(1,1) < datenum(2020,0,0)
        t = localVar(:,1);
        t = t - 0.5/1440;  % subtract 30 seconds from t to put midnight with previous day!
        allCurrentDays = unique(floor(t));
        
    elseif ~isempty(strfind(allFields{i},'NanFlag')) % Nan Flags
        avgTableName = allFields{i}(1:strfind(allFields{i},'NanFlag')-1);
        t = strct.(avgTableName)(:,1);
        t = t - 0.5/1440;  % subtract 30 seconds from t to put midnight with previous day!
        allCurrentDays = unique(floor(t));
        
    elseif ~isempty(strfind(allFields{i},'SpikeFlag')) % Spike Flags
        avgTableName = allFields{i}(1:strfind(allFields{i},'SpikeFlag')-1);
        t = strct.(avgTableName)(:,1);
        t = t - 0.5/1440;  % subtract 30 seconds from t to put midnight with previous day!
        allCurrentDays = unique(floor(t));
        
    elseif strcmp(allFields{i},'rotatedSonic')  % for rotated sonic, use H field
        t = strct.H(:,1);
        t = t - 0.5/1440;  % subtract 30 seconds from t to put midnight with previous day!
        allCurrentDays = unique(floor(t));
    end
    
    % delete days that fall outside fill dates
    for j = 1:numel(allCurrentDays)
        
        if find(fillDates == allCurrentDays(j))
            continue
        else
            localVar(floor(t) == allCurrentDays(j),:) = [];
            t(floor(t) == allCurrentDays(j),:) = [];
        end
    end
    
    % add dates that do not exist in the local variable
    fillMat = [];
    for j = 1:numel(fillDates)
        
        if find(allCurrentDays == fillDates(j))
            continue
        else
            % create missing date with time stamps
            tMat = fillDates(j) + (avgPer:avgPer:1440)./1440;
            nanMat = nan(numRows,numCols);
            nanMat(:,1) = tMat;
            
            % vertically concatnate missing dats
            fillMat = [fillMat; nanMat];
        end
    end
    
    % concatnate nanMat to end of localVar and sort
    localVar = [localVar; fillMat];
    t = [t; fillMat(:,1)];
    [t, order] = sort(t);
    localVar = localVar(order,:);
    
    % delete duplicate rows
    t = round(t*1000);
    [~, index, ~] = unique(t);
    localVar = localVar(index,:);
    
    % store local variable in filledStruct
    filledStrct.(allFields{i}) = localVar;
    
end
end