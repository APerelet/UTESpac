%Will allow multiple files per day

function [data, dataInfo] = loadData(dataFiles,currentDateNumber,totalNumberofDates, info,tableNames)
fprintf('\nEvaluating Date %g of %g in folder %s',currentDateNumber,totalNumberofDates,info.siteFolder);

% create dataInfo cell in base workspace
dataInfo = cell(4,size(dataFiles,2));
data = cell(1,size(dataFiles,2));

% iterate through table files
for ii = 1:size(dataFiles,2)
    % find table name
    tableName = tableNames{ii};

    % find table in siteInfo.m
    tableNumber = strcmpi(tableName,info.tableNames);

    % find table sample frequency
    sampleFrequency = info.tableScanFrequency(tableNumber);

    % find number of expected table columns
    expectedTableColumns = info.tableNumberOfColumns(tableNumber);

    
    % load file into temp
    first = 1;
    clearvars temp tmp
    for jj=1:size(dataFiles, 3)
        if ~isdir(dataFiles{1, ii, jj})
            if first
                % display info
                fprintf('\nloading %s.\nDefined in siteInfo.m\n\tExpected Scan Frequency = %0.02g Hz\n\tNumber of .csv Columns = %0.02g \n\t',dataFiles{ii},sampleFrequency,expectedTableColumns)
                temp = load(dataFiles{1, ii, jj});
                first = 0;
            else
                tmp = load(dataFiles{1, ii, jj});
                temp = [temp; tmp];
            end
        end
    end

    if exist('temp', 'var')
        % check to ensure number of columns is correct, pad or delete if needed
        if size(temp,2) < expectedTableColumns
            warning('number of columns found in %s = %0.02g.  number expected = %0.02g.  Table will be padded with NaNs!',tableName,size(temp,2),expectedTableColumns)
            padding = nan(size(temp,1),expectedTableColumns-size(temp,2));
            temp = horzcat(temp,padding);
        elseif size(temp,2) > expectedTableColumns
            warning('number of columns found in %s = %0.02g.  number expected = %0.02g.  Extra columns will be deleted!',tableName,size(temp,2),expectedTableColumns)
            temp(:,end - (size(temp,2)-expectedTableColumns)+1:end) = [];
        end

        % eliminate round off error in second calculation
        temp(:,4) = round(temp(:,4).*100)./100;

        % sort data in ascending order and ignore duplicates
        dateNumbers = round(temp(:,2).*1e8 + temp(:,3).*10000+temp(:,4)*100);
        [dateNumbers, ia2, ~] = unique(dateNumbers,'R2012a');
        table = temp(ia2,:);

        % find cut off and eliminate all rows beyond 2 days
        cutOff = round((dateNumbers(1)+2*1e8)/100)*100;
        table(dateNumbers>cutOff,:) = [];

        % find table size
        tableSize = size(table);

        % replace -7999 with NaN
        table(table == -7999) = NaN;
        numNans = sum(sum(isnan(table)))/(tableSize(1)*tableSize(2));

        % create empty table with complete time steps (no missed scans!)
        beginSerialDay = datenum(table(2,1),1,0) + table(2,2);
        endSerialDay = datenum(table(end-1,1),1,0) + table(end-1,2);
% % %         if beginSerialDay~=endSerialDay
% % %             error([char(13), 'DEBUG:::Files for day starting on: ', datestr(beginSerialDay),' spans more than 1 day, please check to make sure data file does not contain data going in to next day.', ...
% % %                 char(13), 'This error is for debuging purposes, comment out if datafiles actually span more than 1 day', char(13)]);
% % %         end
        filledTable = completeTableCreate(beginSerialDay,endSerialDay,sampleFrequency,expectedTableColumns);

        % find percent available data
        dataPercent = tableSize(1)/size(filledTable,1);

        % merge table into filledTable - round serial date number to 8 decimals
        [~, ia, ib] = intersect(floor(campbellDate2SerialDate(table(:,1:4))*1e8),floor(campbellDate2SerialDate(filledTable(:,1:4))*1e8),'stable');
        filledTable(ib,2:end) = table(ia,2:end);

        % store filled table in data cell
        data{ii} = filledTable;

        % put table information into dataInfo cell
        localTableInfo{1,1} = sprintf('table %g (%s):',ii,tableName);
        localTableInfo{2,1} = sprintf('\t percent NaNs: %g%s',numNans*100,'%');
        localTableInfo{3,1} = sprintf('\t table scan frequency: %g Hz',sampleFrequency);
        localTableInfo{4,1} = sprintf('\t percent data available: %g%s',dataPercent*100,'%');

        % display results
        for j = 1:length(localTableInfo)
            display(localTableInfo{j,1});
        end

        % assign localTableInfo to dataInfo
        dataInfo{1,1} = 'Run Options';
        dataInfo{2,1} = ['Run on: ',datestr(now,'dd_mmm_yyyy'),' - UTESpac: ',info.UTESpacVersion];
        if strcmp(info.detrendingFormat,'linear')
            dataInfo{3,1} = 'Linear Detrending';
        else
            dataInfo{3,1} = 'Constant Detrending';
        end
        if strcmp(info.PF.globalCalculation,'global')
            dataInfo{4,1} = ['Planar Fit Max/Min Wind: [',num2str(info.PF.globalCalcMinWind),' ',num2str(info.PF.globalCalcMaxWind),']'];
        end
        dataInfo(5:8,ii) = localTableInfo; clear localTableInfo
        fprintf('table %g loaded successfully!\n',ii)
    end
end

% find max number of days spanned by any data table. Possible number of days is 0, 1, or 2
numDays = nan(size(data));
startSerialDay = nan(size(data));
endSerialDay = nan(size(data));

for ii = 1:numel(data)
    
    if isempty(data{ii})
        numDays(ii) = 0;
    else
        startSerialDay(ii) = datenum(data{ii}(1,1),1,0) + data{ii}(1,2);
        endSerialDay(ii) = datenum(data{ii}(end-1,1),1,0) + data{ii}(end-1,2);
        numDays(ii) = endSerialDay(ii) - startSerialDay(ii)+1;
    end
end
[maxNumDaysForAnyTable, index] = max(numDays);

% iterate through all data tables, create non-existent tables and complete incomplete tables
for ii = 1:numel(data)
    
    % find table name
    tableName = tableNames{ii};
        
    % find table in siteInfo.m
    tableNumber = strcmp(tableName,info.tableNames);
        
    if numDays(ii) == 0  % if the table does not exist
        sampleFrequency = info.tableScanFrequency(tableNumber);
        expectedTableColumns = info.tableNumberOfColumns(tableNumber);
        data{ii} = completeTableCreate(startSerialDay(index),endSerialDay(index),sampleFrequency,expectedTableColumns);
    
    elseif numDays(ii) < maxNumDaysForAnyTable % if the table is incomplete 
        if startSerialDay(ii) > startSerialDay(index) % first day is missing
            
            % create filler data for first day
            sampleFrequency = info.tableScanFrequency(tableNumber);
            fillerData = completeTableCreate(startSerialDay(index),startSerialDay(index),sampleFrequency,size(data{ii},2));
                        
            % concatnate filler date with full data
            data{ii} = vertcat(fillerData,data{ii});
            
        elseif endSerialDay(ii) < endSerialDay(index) % second day is missing
            % create filler data for second day
            sampleFrequency = info.tableScanFrequency(tableNumber);
            fillerData = completeTableCreate(endSerialDay(index),endSerialDay(index),sampleFrequency,size(data{ii},2));
                        
            % concatnate filler date with full data
            data{ii} = vertcat(data{ii},fillerData);
        end
    end
end

assignin('base','dataInfo',dataInfo);

end