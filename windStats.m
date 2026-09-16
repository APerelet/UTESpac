function output = windStats(output,sensorInfo,tableNames,info)

% find number sonic info
if isfield(sensorInfo,'u')
    numSonics = size(sensorInfo.u,1);
    
    windEnvelope = info.windDirectionTest.envelopeSize;
    
    % create warning field within output to store warning messages
    output.warnings = cell(1);
    
    for ii = 1:numSonics
            tble = sensorInfo.u(ii,1);
            bearing = sensorInfo.u(ii,4);
            towerBearing = info.tower(ii);
            sonHeight = sensorInfo.u(ii,3);
            uCol = sensorInfo.u(sensorInfo.u(:,3)==sonHeight,2);
            vCol = sensorInfo.v(sensorInfo.v(:,3)==sonHeight,2);
            tableName = tableNames{tble};
            % check sonic manufacturere
            if sensorInfo.u(ii,5)
                u = output.(tableName)(:,uCol);
                v = output.(tableName)(:,vCol);
            else % for RMYoung swap u and v and multiply v by -1
                u = output.(tableName)(:,vCol);
                v = output.(tableName)(:,uCol)*-1;
            end
            t = output.(tableName)(:,1);
            
            % find direction and speed
            dir = mod(atan2(-1*v,u)*180/pi+bearing,360);
            spd = (u.^2+v.^2).^0.5;
            flag = zeros(size(dir));
            
            % flag mean wind speed that falls within envelope
            if towerBearing+windEnvelope > 360
                minAngle = towerBearing - windEnvelope;
                maxAngle = towerBearing-360+windEnvelope;
                flag(dir>minAngle) = 1;
                flag(dir<maxAngle) = 1;
            elseif towerBearing-windEnvelope < 0
                minAngle = 360+towerBearing - windEnvelope;
                maxAngle = towerBearing+windEnvelope;
                flag(dir>minAngle) = 1;
                flag(dir<maxAngle) = 1;
            else
                minAngle = towerBearing - windEnvelope;
                maxAngle = towerBearing + windEnvelope;
                temp1 = zeros(size(dir));
                temp2 = zeros(size(dir));
                temp1(dir>minAngle) = 1;
                temp2(dir<maxAngle) = 1;
                flag((temp1+temp2) == 2) = 1;
            end
            output.spdAndDir(1:length(spd),1) = t;
            output.spdAndDir(1:length(spd),ii*3-1:ii*3+1) = [dir spd flag];
            output.spdAndDirHeader{1} = 'time';
            output.spdAndDirHeader{ii*3-1} = strcat(num2str(sonHeight), 'm Direction [deg]');
            output.spdAndDirHeader{ii*3} = strcat(num2str(sonHeight), 'm Speed [m/s]');
            output.spdAndDirHeader{ii*3+1} = strcat(num2str(sonHeight), 'm Dir_Flag_', num2str(minAngle), '<Dir<', num2str(maxAngle), '[-]');
    end
end
