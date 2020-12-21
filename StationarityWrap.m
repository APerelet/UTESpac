function [output] = StationarityWrap(data, rotatedSonicData, info, output, sensorInfo, tableNames)

StationarityHeader = {'Timestamp'; []};
velComp = {'u', 'v', 'w'};

if isfield(sensorInfo, 'u')
    numSonics = size(sensorInfo.u, 1);
    
    % Iterate through sonics
    for ii=1:numSonics

        %Calculate number of points per averaging period - assume all sonics on
        %same table...will need to fix this later
        TableNum = sensorInfo.u(ii, 1);
        sonHeight = sensorInfo.u(ii, 3); 
        tmp = regexp(info.tableNames, tableNames{TableNum});
        [~, siteInfoTableNum] = max(~cellfun(@isempty, tmp));
        pnts = info.tableScanFrequency(:, siteInfoTableNum)*60*info.avgPer;

        M = pnts/5;   % Set as option later

        % Create average time vector
        time = data{TableNum}(pnts:pnts:end, 1);

        cntr1 = 1;
        cntr2 = pnts;
        for qq=1:length(data{TableNum})/pnts

            if ii==1
                colIndStart = 1;
            else
                colIndStart = colInd;
                
            end
            colInd = colIndStart;
            
            % finewire temperature
            if isfield(sensorInfo, 'fw')

                %Temperature
                X = data{TableNum}(cntr1:cntr2, sensorInfo.fw(ii, 2));
                Y = data{TableNum}(cntr1:cntr2, sensorInfo.fw(ii, 2));
                [Flag(qq, colInd), ~] = ...
                    Stationarity(X, Y, M);
                
                % Add information to header
                if qq==1
                    StationarityHeader = [StationarityHeader, {'T_fw Flag'; num2str(sensorInfo.fw(ii, 3))}];
                end
                colInd = colInd + 1;

                %Sensible Heat Flux
                X = rotatedSonicData(cntr1:cntr2, 3*ii);
                Y = data{TableNum}(cntr1:cntr2, sensorInfo.fw(ii, 2));
                [Flag(qq, colInd), ~] = ...
                    Stationarity(X, Y, M);

                % Add information to header
                if qq==1
                    StationarityHeader = [StationarityHeader, {'w''t_fw'' Flag'; num2str(sensorInfo.fw(ii, 3))}];
                end

                colInd = colInd + 1;
            end


            % Sonic temperature
            if isfield(sensorInfo, 'Tson')
                % Temperature
                X = data{TableNum}(cntr1:cntr2, sensorInfo.Tson(ii, 2));
                Y = data{TableNum}(cntr1:cntr2, sensorInfo.Tson(ii, 2));
                [Flag(qq, colInd), ~] = ...
                    Stationarity(X, Y, M);
                
                % Add information to header
                if qq==1
                    StationarityHeader = [StationarityHeader, {'Tson Flag'; num2str(sensorInfo.Tson(ii, 3))}];
                end
                colInd = colInd + 1;

                % Sensible Heat Flux
                X = rotatedSonicData(cntr1:cntr2, 3*ii);
                Y = data{TableNum}(cntr1:cntr2, sensorInfo.Tson(ii, 2));
                [Flag(qq, colInd), ~] = ...
                    Stationarity(X, Y, M);

                % Add information to header
                if qq==1
                    StationarityHeader = [StationarityHeader, {'w''t_son'' Flag'; num2str(sensorInfo.Tson(ii, 3))}];                
                end
                colInd = colInd + 1;    
            end

            % Moisture and Latent Heat Flux
            if isfield(sensorInfo, 'irgaH2O')
                checkHeight = find(sensorInfo.irgaH2O(:, 3)==sonHeight);
                
                % Moisture
                X = data{TableNum}(cntr1:cntr2, sensorInfo.irgaH2O(checkHeight, 2));
                Y = data{TableNum}(cntr1:cntr2, sensorInfo.irgaH2O(checkHeight, 2));
                [Flag(qq, colInd), ~] = ...
                    Stationarity(X, Y, M);
                
                % Add information to header
                if qq==1
                    StationarityHeader = [StationarityHeader, {'H2O Flag'; num2str(sensorInfo.irgaH2O(checkHeight, 3))}];
                end
                
                colInd = colInd + 1;
                
                % Latent Heat Flux
                X = rotatedSonicData(cntr1:cntr2, 3*ii);
                Y = data{TableNum}(cntr1:cntr2, sensorInfo.irgaH2O(checkHeight, 2));
                [Flag(qq, colInd), ~] = ...
                    Stationarity(X, Y, M);

                % Add information to header
                if qq==1
                    StationarityHeader = [StationarityHeader, {'w''q'' Flag'; num2str(sensorInfo.irgaH2O(checkHeight, 3))}];                
                end
            end

            % update counters
            cntr1 = cntr1 + pnts;
            cntr2 = cntr2 + pnts; 
        end
    end

    output.Stationarity = [time, Flag];
    output.StationarityHeader = StationarityHeader;
else
    % No sonics were found
    warning([char(10), '----------------------------------', char(10), char(13), 'No Sonics Found...', char(10), char(13), '----------------------------------']);
end
