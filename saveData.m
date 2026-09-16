function output = saveData(info,output,dataInfo,tableNames, rawFlux, template, fileDate)
% saveData saves output data in the site folder withinn the root folder
saveDir = [info.rootFolder, filesep, info.siteFolder, info.foldStruct];
fprintf('\nSaving data...\n');

outputDir = ['output', num2str(info.avgPer)];

if ~exist([saveDir, filesep, outputDir], 'dir')
    fprintf(['\n\tOutput directory for', num2str(info.avgPer), ' Min average does not exist:\nCreating output dicrectory.\n']);
    mkdir([saveDir, filesep, outputDir]);
end
% include dataInfo and headers in output folder
output.dataInfo = dataInfo;
output.tableNames = tableNames;

% sort structure
output = orderfields(output);

% find correct filename and save output
if strcmp(info.PF.globalCalculation,'global')
    PFtype = 'GPF_';
else
    PFtype = 'LPF_';
end
if strcmp(info.detrendingFormat,'linear')
    detrendType = 'LinDet_';
elseif strcmp(info.detrendingFormat,'constant')
    detrendType = 'ConstDet_';
elseif strcmp(info.detrendingFormat,'wavelet')
    detrendType = 'MRADet_';
end
%% save .mat structure
fileName = strcat(info.siteFolder(5:end),'_',num2str(info.avgPer),'minAvg_',PFtype,detrendType,fileDate,'.mat');
save(strcat(saveDir,filesep,outputDir, filesep, fileName),'output', '-v7.3');
%% save .csv file
if info.saveCSV
    csvSave(template,output,info)
end
%% save netCDF
if info.saveNetCDF

    netCDFtDir = 'outputNetCDF';

    if ~exist([saveDir, filesep, netCDFtDir], 'dir')
        fprintf('\n\tOutput directory does not exist: Creating output dicrectory.\n');
        mkdir([saveDir, filesep, netCDFtDir]);
    end
    
    fileName = strcat(info.siteFolder(5:end),'_',num2str(info.avgPer),'minAvg_',PFtype,detrendType,fileDate,'.nc');
    
    ncFullFileName = strcat(saveDir,filesep,netCDFtDir,filesep,fileName);
    
    % Delete file if it exists
    delete(ncFullFileName);

    fields = fieldnames(output);

    %Extract Headers
    headerFlag = cellfun(@(x) contains(x, 'Header'), fields);   % Find headers
    tableNameFlag = cellfun(@(x) contains(x, 'tableNames'), fields);    % Find table name field
    tableFlag = cellfun(@(x) contains(x, 'Flag'), fields);              % Find Flag fields
    warningFlag = cellfun(@(x) contains(x, 'warning'), fields);         % Find warning field
    StructFuncFlag = cellfun(@(x) contains(x, 'StructFunc'), fields);         % Find warning field
    
    % Header Fields
    fieldHeaders = fields(headerFlag);
    % Data Fields
    fieldData = fields(~headerFlag & ~tableNameFlag & ~tableFlag & ~warningFlag & ~StructFuncFlag);
    
    % Iterate through output structure Data fields
    for ii=1:length(fieldData)
        % Check if header exists
        CurrentHeaderFlag = cellfun(@(x) startsWith(x, fieldData{ii}), fieldHeaders);
    
        % Check if field is from raw data table and correct header
        rawTableFlag = cellfun(@(x) strcmp(fieldData{ii}, x), output.tableNames);
    
        if any(rawTableFlag)
            tmpHeader = {'time'};
            for qq=2:length(output.([output.tableNames{rawTableFlag}, 'Header']))
    
                tmpHeader{qq} = [num2str(output.([output.tableNames{rawTableFlag}, 'Header']){2, qq}), 'm ', ...
                    output.([output.tableNames{rawTableFlag}, 'Header']){1, qq}, ' [-]'];
            end
            output.([output.tableNames{rawTableFlag}, 'Header']) = tmpHeader;
        end
    
        if any(CurrentHeaderFlag)
            % Check if header Exists
            for jj = 1:length(output.(fieldHeaders{CurrentHeaderFlag}))
    
                % Check Data columns and save netCDF entry
                if strcmp(output.(fieldHeaders{CurrentHeaderFlag}){jj}, 'time')
                    % Check for time variable
                    if exist(ncFullFileName, "file")
                        tmp = ncinfo(ncFullFileName);
                        timeFlag = any(cellfun(@(x) contains(x, 'time'), {tmp.Variables.Name}));
                    else
                        timeFlag = 0;
                    end
                    varName = 'time';
                    varUnits = 'Time';
                    varData = string(datestr(output.(fieldData{ii})(:, jj), 'yyyy mm dd HH:MM:SS.FFF'));
                    varDim = size(varData);
    
                    % Create NC entry for time
                    if ~timeFlag
                        nccreate(ncFullFileName, varName, "Dimensions", {'x', varDim(1), 'y', varDim(2)}, "Datatype","string", "FillValue", "NaN", "Format","netcdf4");
                        ncwrite(ncFullFileName, varName, varData);
                        ncwriteatt(ncFullFileName, varName, 'Units', varUnits, 'Datatype','string');
                    end
                
                else
                    % Parse height, variable name, and units from data header
                    varInfo = regexp(output.(fieldHeaders{CurrentHeaderFlag}){jj}, '(?<Height>(\d*.?\d*|N\/A)m)\s*(?<VarName>[^\[\]]+)\s?\[(?<Units>.*)\]', 'names');
        
                    varName = [varInfo.VarName, '_', varInfo.Height(1:end-1)];
                    varUnits = varInfo.Units;
                    varHeight = varInfo.Height;
                    varData = output.(fieldData{ii})(:, jj);
                    varDim = size(varData);
        
                    % Create NC Entry for data column
    
                    % Check if varName Exists
                    tmp = ncinfo(ncFullFileName);
                    if any(cellfun(@(x) contains(varName, x), {tmp.Variables.Name}))
                        warning(['Found Duplicate Variable with name ', varName, ' in ', fieldHeaders{CurrentHeaderFlag}, '. Skipping duplicate instance. Fix header if this variable is important']);
                    else
                        nccreate(ncFullFileName, varName, "Dimensions", {'x', varDim(1), 'y', varDim(2)}, "Datatype","double", "FillValue", NaN, "Format","netcdf4");
                        ncwrite(ncFullFileName, varName, varData);
                        ncwriteatt(ncFullFileName, varName, 'Units', varUnits, 'Datatype','string');
                        ncwriteatt(ncFullFileName, varName, 'Height', varHeight, 'Datatype','string');
                    end
                end
            end
        else
            % If no header exists do not save data
            warning(['No header found for field ', fieldData{ii}, '. Skipping field. If important please include header!']);
        end
    end
     
end
%% save raw fluxes
if info.saveRawConditionedData
    
    rawDir = 'outputRAW';
    
    if ~exist([saveDir, filesep, rawDir], 'dir')
        fprintf('\n\tOutput directory for RAW does not exist:\nCreating Raw output dicrectory.\n');
        mkdir([saveDir, filesep, rawDir]);
    end
    
    fileName = strcat(info.siteFolder(5:end),'_raw_',PFtype,detrendType,fileDate);
    save(strcat(saveDir,filesep,rawDir,filesep,fileName),'rawFlux', '-v7.3');
end
end
