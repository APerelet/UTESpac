function [output] = StationarityWrap(data, rotatedSonicData, info, output, sensorInfo, tableNames)

StationarityHeader = {'Timestamp'; []};
velComp = {'u', 'v', 'w'};


%Calculate number of points per averaging period - assume all sonics on
%same table...will need to fix this later
TableNum = sensorInfo.u(1, 1);
tmp = regexp(info.tableNames, tableNames{TableNum});
[~, siteInfoTableNum] = max(~cellfun(@isempty, tmp));
pnts = info.tableScanFrequency(:, siteInfoTableNum)*60*info.avgPer;
    
M = pnts/5;   % Set as option later - since length of data is 1 day, set stationary test to 5 minute chunks


% Create average time vector
time = data{TableNum}(pnts:pnts:end, 1);

cntr1 = 1;
cntr2 = pnts;
for qq=1:length(data{TableNum})/pnts
    
    colInd = 1;
    % finewire temperature
    if isfield(sensorInfo, 'fw')
        for ii=1:size(sensorInfo.fw, 1)

% % %             %Temperature
% % %             X = data{TableNum}(cntr1:cntr2, sensorInfo.fw(ii, 2));
% % %             Y = data{TableNum}(cntr1:cntr2, sensorInfo.fw(ii, 2));
% % %             [Flag(qq, colInd), ~] = ...
% % %                 Stationarity(X, Y, M);
% % %             
% % %             % Add information to header
% % %             StationarityHeader = [StationarityHeader, {'T_fw Flag'; num2str(sensorInfo.fw(ii, 3))}];
% % % 
% % %             colInd = colInd + 1;

            %Sensible Heat Flux
            X = rotatedSonicData(cntr1:cntr2, 3);
            Y = data{TableNum}(cntr1:cntr2, sensorInfo.fw(ii, 2));
            [Flag(qq, colInd), ~] = ...
                Stationarity(X, Y, M);
            
            % Add information to header
            if qq==1
                StationarityHeader = [StationarityHeader, {'w''t_fw'' Flag'; num2str(sensorInfo.fw(ii, 3))}];
            end

            colInd = colInd + 1;
        end
    end


    % Sonic temperature
    if isfield(sensorInfo, 'Tson')
        for ii=1:size(sensorInfo.Tson, 1)
% % %             % Temperature
% % %             X = data{TableNum}(cntr1:cntr2, sensorInfo.Tson(ii, 2));
% % %             Y = data{TableNum}(cntr1:cntr2, sensorInfo.Tson(ii, 2));
% % %             [Flag(qq, colInd), ~] = ...
% % %                 Stationarity(X, Y, M);
% % %             
% % %             % Add information to header
% % %             StationarityHeader = [StationarityHeader, {'Tson Flag'; num2str(sensorInfo.Tson(ii, 3))}];
% % % 
% % %             colInd = colInd + 1;

            % Sensible Heat Flux
            X = rotatedSonicData(cntr1:cntr2, 3);
            Y = data{TableNum}(cntr1:cntr2, sensorInfo.Tson(ii, 2));
            [Flag(qq, colInd), ~] = ...
                Stationarity(X, Y, M);
            
            % Add information to header
            if qq==1
                StationarityHeader = [StationarityHeader, {'w''t_son'' Stationary Flag'; num2str(sensorInfo.Tson(ii, 3))}];                
            end
            colInd = colInd + 1;    
        end
    end

    % Moisture and Latent Heat Flux
    
    % update counters
    cntr1 = cntr1 + pnts;
    cntr2 = cntr2 + pnts; 
end

output.Stationarity = [time, Flag];
output.StationarityHeader = StationarityHeader;
