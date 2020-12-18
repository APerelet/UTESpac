%Allows for multiple files per day

function [headersCell, dataFiles, tableNames, info] = findFiles(info)

% find possible sites of interest
sitesStruct = dir(strcat(info.rootFolder, filesep,'site*'));

% display possible sites to command window
for ii = 1:numel(sitesStruct)
    display(sprintf('%g. %s',ii,sitesStruct(ii).name));
end

% ask user to select appropriate site
site = sitesStruct(input('Please indicate site number of interest: ')).name; clc;
info.siteFolder = site;

% store site folder
siteFolder = strcat(info.rootFolder,filesep,site, info.foldStruct);

% run siteInfo script
run(strcat(siteFolder,filesep,'siteInfo.m'));

% find headers for selected site
headers = dir(strcat(siteFolder,filesep,'*header.*'));



if isempty(headers)
    error('Unable to Find Table Headers!  Check format of header name.')
end

% iterate through headers and store .csv files in csvFileStruct
for ii = 1:length(headers)
    currentHeaderFileName = headers(ii).name;
    
    % Ignore Tables
    if isfield(info, 'TableIgnore')
        if any(cell2mat(cellfun(@(x) ~isempty(strfind(currentHeaderFileName, x)), info.TableIgnore, 'UniformOutput', 0)))
            fprintf(['\n---------------------\nSkipping files associated with: ',...
                currentHeaderFileName, '\n---------------------\n']);
            ignoreFlag(ii) = 1;
            continue
        else
            ignoreFlag(ii) = 0;
        end
    else
        ignoreFlag(ii) = 0;
    end

    
    
    try
        % find current header array
        headerFile = fopen(strcat(siteFolder,filesep,currentHeaderFileName));
        headerLine = fgetl(headerFile);
        currentHeaderArray =  strsplit(headerLine, ',');
               
        % store measurement height on second row of currentHeaderArray
        for jj = 1:length(currentHeaderArray)
            
            % find digits in column title
            digitLocations = find(isstrprop(currentHeaderArray{1,jj},'digit')==1);
            
            % find decimals in column title
            decimalLocations = strfind(currentHeaderArray{1,jj},'.');
            
            % concatnate digits and decimals, sort in ascending order
            digitAndDecimalLocations = sort(horzcat(digitLocations, decimalLocations));
            
            % if digitsAndDecimals is empty, no height was found 
            if isempty(digitAndDecimalLocations)
                currentHeaderArray{2,jj} = [];
                
                % if digitsAndDecimals has 1 value, this must be the height 
            elseif numel(digitAndDecimalLocations) == 1
                heightLocations = digitAndDecimalLocations;
                currentHeaderArray{2,jj} = str2double(currentHeaderArray{1,jj}(heightLocations));
                
            else % if digitsAndDecimals has multiple values, use diff to find consecutive locations, 
                 % last chunk of consecutive locations must be the height locations
                heightBeginMarker = find(diff(digitAndDecimalLocations) > 1,1,'last');
                if isempty(heightBeginMarker)
                    heightBeginMarker = 0;
                end
                heightLocations = digitAndDecimalLocations(heightBeginMarker+1:end);
                currentHeaderArray{2,jj} = str2double(currentHeaderArray{1,jj}(heightLocations));
            end
        end
        
        
        % store headers in cell 
        headersCell{ii} = currentHeaderArray; clear currentHeaderArray
        
    catch err
        error('Unable to load %s \@ ln %g.  Check format! \n%s',currentHeaderFileName,err.message);
    end
    
    % find length of table name
    tableNameLength = strfind(currentHeaderFileName,'_header')-1;
    
    % remove '_header' to find table name
    tableName = currentHeaderFileName(1:tableNameLength);
    
    % store table names in cell
    tableNames{ii} = tableName;
    
    % find all files associated with the tableName
    tableFiles = dir(strcat(siteFolder,filesep,strcat('*',tableName,'*')));

    % eliminate the header file from the list
    tableFiles(~cellfun(@isempty,strfind({tableFiles(:).name},'header'))) = [];
    
    if isempty(tableFiles)
       error(['Could not find table ', tableName, char(10), ...
           char(9), 'Please confirm header and data files have same name', char(10)]);
    end
    
    %In case of similar table names, omit wrong table names
    clearvars nameFlag
    for kk = 1:length(tableFiles)
        tableFile = tableFiles(kk).name;

        %Parse filename to get date
        fileForm = regexp(tableFile, info.FileForm, 'names');
        FileFormFlag(kk) = 0;
        if isempty(fileForm)
            FileFormFlag(kk) = 1;
            info.FileForm1 = [info.FileForm(1:end-4), '_\d+.dat'];
            fileForm = regexp(tableFile, info.FileForm1, 'names');
            if isempty(fileForm)
                error(['Invalid Fileform for file: ', char(10), ...
                    char(9), tableFile, char(10), ...
                    'Please insure only valid table files and header files exist in directory: ', char(10),...
                    char(9), tableFiles(kk).folder, char(10),...
                    'If file is valid please make sure the naming satisfies the regular expression defined in info.FileForm', char(10)]);
            end
        end
        
        nameCheck = strcmp(fileForm.TableName, tableName);
        if nameCheck
            nameFlag(kk) = true;
        else
            nameFlag(kk) = false;
        end
    end
    tableFiles = tableFiles(nameFlag);
    
    
    % sort and store csv files in csvFilesCell
    csvFiles = cell(length(tableFiles),1);  %preallocate for speed
    csvDates = NaN(size(csvFiles));
    for jj = 1:length(tableFiles)
        tableFile = tableFiles(jj).name;

        %Parse filename to get date
        if FileFormFlag(jj)
            fileForm = regexp(tableFile, info.FileForm1, 'names');
        else
            fileForm = regexp(tableFile, info.FileForm, 'names');
        end
        
        %Check if Day Hour Minutes exist
        testString = {'Day', 'Hour', 'Minute', 'Second'};
        fieldnames = fields(fileForm);
        contentCheck = cell2mat(cellfun(@(x) strcmp(x, testString), fieldnames, 'UniformOutput', false));
        
        %Each filename must contain a Day
        if ~any(contentCheck(:, 1))
            error('File format needs to have at least Year, Month, and Day. Check input for fileForm and change file names if necessary');
        end
        
        %If no Hour in filename, set to 0
        if ~any(contentCheck(:, 2))
            fileForm.Hour = '0';
        end
        %If no minutes in filename, set to 0
        if ~any(contentCheck(:, 3))
            fileForm.Minute = '0';
        end

        %If no minutes in filename, set to 0
        if ~any(contentCheck(:, 4))
            fileForm.Second = '0';
        end
        
        
        % place tableFile name in csvFiles cell
        csvFiles{jj} = tableFile;
        
        % parse date string to create date number
        csvDates(jj) = datenum(str2double(fileForm.Year), str2double(fileForm.Month), str2double(fileForm.Day), str2double(fileForm.Hour), str2double(fileForm.Minute), str2double(fileForm.Second));
    end
    
    %Find Unique Days
    [a, ind] = unique(floor(csvDates));
    ind(end+1) = length(csvDates);
    
    % place csv files in cell matrix where the row is determined by the date
    if ii == 1
        dateBegin = floor(min(csvDates)-20); % all tables must begin and end within 20 days of first and last dates of table1
        dateEnd = floor(max(csvDates)+20);
% % %         if length(ind)==2
% % %             dataFiles = cell(1,length(headers));
% % %         else
            dataFiles = cell(dateEnd-dateBegin,length(headers), round(max(diff(ind))/10)*10);
% % %         end
    end
    
    % Check if there are files with outlier date stamps
    dateCheck = find(or(a<dateBegin, a>dateEnd));
        if ~isempty(dateCheck)
           error(['The following Files do not fall outside expected time:', char(10),...
               csvFiles{ind(dateCheck)}, char(10),...
               'Please Delete or fix the dates']);
        end
    
    if length(ind)==2
        dataFiles{1, ii, 1} = csvFiles{1};
    else
        for jj = 1:(length(ind)-1)
            for qq=1:length(ind(jj):ind(jj+1)-1)
                dataFiles{a(jj)-dateBegin,ii, qq} = csvFiles{ind(jj)+qq-1};
            end
        end
    end
end

%Remove Ignored tables
headersCell = headersCell(~ignoreFlag);
dataFiles = dataFiles(:, ~ignoreFlag, :);
tableNames = tableNames(~ignoreFlag);


% check cell matrix for empty rows and delete them
dataFilesFlag = any((sum(cellfun(@isempty,dataFiles),3)==size(dataFiles, 3))==0, 2);

dataFiles = dataFiles(dataFilesFlag, :, :);

% ask user to select which dates to calculate
displayCell = cell(size(dataFiles,1),size(dataFiles,2)+1);
displayCell(1:end,1) = num2cell(1:size(dataFiles,1))';
displayCell(1:end,2:end) = dataFiles(:, :, 1);

fprintf('\nDisplaying all files.\nNote: only displaying first file found for each day.');
showcell(displayCell);

SOI = input('Plese input dates of interest. e.g. [1 3 4:7] or ''0'' for all dates: ');

if SOI > 0
    dataFiles = dataFiles(SOI,:, :);
end
for ii = 1:size(dataFiles,2)
    for jj = 1:size(dataFiles,1)
        for kk=1:size(dataFiles, 3)
            dataFiles{jj,ii, kk} = strcat(siteFolder,filesep,dataFiles{jj,ii, kk});
        end
    end
end
clc
end


