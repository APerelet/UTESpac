function [output] = StationarityWrap(data, rotatedSonicData, info, output, sensorInfo, tableNames)

StationarityHeader = {'time'};
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

            if qq==1
                if ii==1
                    colIndStart = 1;
                else
                    colIndStart = colInd;

                end
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
                    StationarityHeader = [StationarityHeader, {[num2str(sensorInfo.fw(ii, 3)), 'm Tfw_Stationarity_Flag [-]']}];
                end
                colInd = colInd + 1;

                %Sensible Heat Flux
                X = rotatedSonicData(cntr1:cntr2, 3*ii);
                Y = data{TableNum}(cntr1:cntr2, sensorInfo.fw(ii, 2));
                [Flag(qq, colInd), ~] = ...
                    Stationarity(X, Y, M);

                % Add information to header
                if qq==1
                    StationarityHeader = [StationarityHeader, {[num2str(sensorInfo.fw(ii, 3)), 'm KinHeatFlux_wPF''Tfw''_Stationarity_Flag [-]']}];
                end

                colInd = colInd + 1;
            end


            % Sonic temperature
            if isfield(sensorInfo, 'Tson')
		%Velocity Components
                % Streamwise
                X = data{TableNum}(cntr1:cntr2, sensorInfo.u(ii, 2));
                Y = data{TableNum}(cntr1:cntr2, sensorInfo.u(ii, 2));
                [Flag(qq, colInd), ~] = ...
                    Stationarity(X, Y, M);
                
                % Add information to header
                if qq==1
                    StationarityHeader = [StationarityHeader, {[num2str(sensorInfo.u(ii, 3)), 'm Streamwise_u_Stationarity_Flag [-]']}];
                end
                colInd = colInd + 1;
                
                % Spanwise
                X = data{TableNum}(cntr1:cntr2, sensorInfo.v(ii, 2));
                Y = data{TableNum}(cntr1:cntr2, sensorInfo.v(ii, 2));
                [Flag(qq, colInd), ~] = ...
                    Stationarity(X, Y, M);
                
                % Add information to header
                if qq==1
                    StationarityHeader = [StationarityHeader, {[num2str(sensorInfo.u(ii, 3)), 'm Spanwise_v_Stationarity_Flag [-]']}];
                end
                colInd = colInd + 1;
                
                % Vertical
                X = data{TableNum}(cntr1:cntr2, sensorInfo.w(ii, 2));
                Y = data{TableNum}(cntr1:cntr2, sensorInfo.w(ii, 2));
                [Flag(qq, colInd), ~] = ...
                    Stationarity(X, Y, M);
                
                % Add information to header
                if qq==1
                    StationarityHeader = [StationarityHeader, {[num2str(sensorInfo.u(ii, 3)), 'm Vertical_w_Stationarity_Flag [-]']}];
                end
                colInd = colInd + 1;

                % Temperature
                X = data{TableNum}(cntr1:cntr2, sensorInfo.Tson(ii, 2));
                Y = data{TableNum}(cntr1:cntr2, sensorInfo.Tson(ii, 2));
                [Flag(qq, colInd), ~] = ...
                    Stationarity(X, Y, M);
                
                % Add information to header
                if qq==1
                    StationarityHeader = [StationarityHeader, {[num2str(sensorInfo.u(ii, 3)), 'm TSon_Stationarity_Flag [-]']}];
                end
                colInd = colInd + 1;

                % Sensible Heat Flux
                X = rotatedSonicData(cntr1:cntr2, 3*ii);
                Y = data{TableNum}(cntr1:cntr2, sensorInfo.Tson(ii, 2));
                [Flag(qq, colInd), ~] = ...
                    Stationarity(X, Y, M);

                % Add information to header
                if qq==1
                    StationarityHeader = [StationarityHeader, {[num2str(sensorInfo.Tson(ii, 3)), 'm KinHeatFlux_wPF''TSon''_Stationarity_Flag [-]']}];             
                end
                colInd = colInd + 1;    
            end

            % Moisture and Latent Heat Flux
            if isfield(sensorInfo, 'irgaH2O')
                checkHeight = find(sensorInfo.irgaH2O(:, 3)==sonHeight);
                
                if ~isempty(checkHeight)
                    % Moisture
                    X = data{TableNum}(cntr1:cntr2, sensorInfo.irgaH2O(checkHeight, 2));
                    Y = data{TableNum}(cntr1:cntr2, sensorInfo.irgaH2O(checkHeight, 2));
                    [Flag(qq, colInd), ~] = ...
                        Stationarity(X, Y, M);

                    % Add information to header
                    if qq==1
                        StationarityHeader = [StationarityHeader, {[num2str(sensorInfo.irgaH2O(checkHeight, 3)), 'm Humidity_Stationarity_Flag [-]']}];
                    end

                    colInd = colInd + 1;

                    % Latent Heat Flux
                    X = rotatedSonicData(cntr1:cntr2, 3*ii);
                    Y = data{TableNum}(cntr1:cntr2, sensorInfo.irgaH2O(checkHeight, 2));
                    [Flag(qq, colInd), ~] = ...
                        Stationarity(X, Y, M);

                    % Add information to header
                    if qq==1
                        StationarityHeader = [StationarityHeader, {[num2str(sensorInfo.irgaH2O(checkHeight, 3)), 'm KinHumidityFlux_Stationarity_Flag [-]]']}];             
                    end
                    
                    colInd = colInd + 1;
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
