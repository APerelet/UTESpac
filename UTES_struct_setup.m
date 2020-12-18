%Prepare to calculate 2nd order structure parameters of temperature, humidity, and
%the temperature humidity cross structure parameter

%Can add to calculate further structure parameters here at a later date.

function [output, raw] = UTES_struct_setup(data, info, output, raw, tableNames, sensorInfo)

if isempty(raw)
    warning('Could not calculate Structure parameters. raw data does not exist')
    return
end

fprintf('\nCalculating Structure Parameters\nC_T^2 C_Tq C_q^2\n');

Rd = 287.058;  % [J/K/kg] Gas constant for air


%%%%%%%%%%%%%
%When Finewire does not exist
%Eventually need a better system than this
fields = fieldnames(sensorInfo);
no_fw = ~any(cellfun(@(x) strcmp(x, 'fw'), fields));
no_irga = ~any(cellfun(@(x) strcmp(x, 'irgaH2O'), fields));
%%%%%%%%%%%%%


%Number of Sonics
NumSonic = size(sensorInfo.Tson, 1);
%Number IRGAs
if no_irga
    NumIRGA = 0;
else
    NumIRGA = size(sensorInfo.irgaH2O, 1);
end

MaxInst = max([NumSonic, NumIRGA]);
MinInst = min([NumSonic, NumIRGA]);



%Iterate through number of instruments (heights & Tables)
cntrIRGA = 1;
for qq=1:MaxInst
    %TableName columns in UTESpac and siteInfo do not match....
    %These lines look for the tablename from siteInfo
    TableNum = sensorInfo.Tson(qq, 1);
    tmp = regexp(info.tableNames, tableNames{TableNum});
    [~, siteInfoTableNum] = max(~cellfun(@isempty, tmp));
    
    % measurement frequency
    freq = info.tableScanFrequency(siteInfoTableNum);
    z = sensorInfo.u(qq, 3);

    %Calculate number of points per averaging period
    pnts = info.tableScanFrequency(:, siteInfoTableNum)*60*info.avgPer;

    %Create average time vector
    time = raw.t(pnts:pnts:end);
   
    %Iterate through each aeraging period and calculate strucutre
    %parameters
    cntr1=1;
    cntr2 = pnts;
    for ii=1:length(raw.t)/pnts
        
        % find reference pressure
        zRef = min(sensorInfo.u(:,3));  % zRef is lowest sonic level
        altitude = info.siteElevation;
        PRef = 101325*(1-2.25577*10^-5*(altitude+zRef))^5.25588; % find pRef from elevation: http://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html

        SonHeight = sensorInfo.Tson(qq, 3);
        %Sonic Variables for each height
        sigma_u = nanstd(raw.uPF(cntr1:cntr2, qq)); 
        sigma_v = nanstd(raw.vPF(cntr1:cntr2, qq)); 
        sigma_w = nanstd(raw.wPF(cntr1:cntr2, qq)); 
        U = nanmean(sqrt(raw.uPF(cntr1:cntr2, qq).^2+raw.vPF(cntr1:cntr2, qq).^2));
        if isfield(sensorInfo, 'P')
            if and(size(sensorInfo.P, 1)==1, sensorInfo.P(1)==siteInfoTableNum)
                P = nanmean(data{TableNum}(cntr1:cntr2, sensorInfo.P(1, 2)).*1000); %average Pressure [Pa]
            else
                P = nanmean(data{TableNum}(cntr1:cntr2, sensorInfo.P(siteInfoTableNum, 2)).*1000); %average Pressure [Pa]
            end

            if isnan(P) || abs((P - PRef)/PRef) > 0.05
                P = PRef;
            end
        else
            P = PRef;
        end

        if ~no_fw
            T = raw.fwT(cntr1:cntr2, qq)+273.15; %Temperature [K]
            
            %Calculate temperature structure parameter
            [Dt2(ii, :), Ct2(ii)] = ...
                crossStruct(T, T, freq, z, info.StructParamSep, info.StructParam_rDiff, 'spatial', sigma_u, sigma_v, sigma_w, U); 
            if Ct2(ii)<0
                Ct2(ii) = nan;
            end
        else
            T = nan;
        end
        TSon = raw.sonTs(cntr1:cntr2, qq)+273.15; %Sonic Temperature [k]
        
        
        %Calculate temperature structure parameter with Sonic temperature
        [Dtson2(ii, :), Ctson2(ii)] = ...
            crossStruct(TSon, TSon, freq, z, info.StructParamSep, info.StructParam_rDiff, 'spatial', sigma_u, sigma_v, sigma_w, U); 
        if Ctson2(ii)<0
            Ctson2(ii) = nan;
        end

        
        if cntrIRGA<=MinInst
            %look for matching IRGA for Sonic height
            IRGAHeight = find(sensorInfo.irgaH2O(:, 3)==SonHeight);

            if ~isempty(IRGAHeight)
                %Abslute humidity [kg m^-3]
                rhov = raw.rhov(cntr1:cntr2, IRGAHeight)./1000;  

                %Calculate mean temperature for period
                T_mean = nanmean(T);
                if isnan(T_mean)
                    T_mean = nanmean(TSon);
                end

                %Specific humidity [kg/kg]
                q = (rhov./(P./(Rd*T_mean)));

                %Calculate temperature humidity cross structure parameter
                
                
                %Calculate humidity structure parameter
                [Dq2(ii, :), Cq2(ii)] = ...
                    crossStruct(q, q, freq, z, info.StructParamSep, info.StructParam_rDiff, 'spatial', sigma_u, sigma_v, sigma_w, U);
                if Cq2(ii)<0
                    Cq2(ii) = nan;
                end

                %Calculate temperature humidity correlation
                if ~no_fw
                    [Dtq(ii, :), Ctq(ii)] = ...
                        crossStruct(T, q, freq, z, info.StructParamSep, info.StructParam_rDiff, 'spatial', sigma_u, sigma_v, sigma_w, U);
                    
                    r_tq(ii) = Ctq(ii)./sqrt(Ct2(ii).*Cq2(ii));
                end
                
                %Calculate temperature humidity correlation TSonic
                [Dtsonq(ii, :), Ctsonq(ii)] = ...
                    crossStruct(TSon, q, freq, z, info.StructParamSep, info.StructParam_rDiff, 'spatial', sigma_u, sigma_v, sigma_w, U);
                
                %Calculate temperature humidity correlation
                r_tsonq(ii) = Ctsonq(ii)./sqrt(Ctson2(ii).*Cq2(ii));

                IRGAFlag = 1;
            else
                IRGAFlag = 0;
            end
        else
            IRGAFlag = 0;
        end
        %update counters
        cntr1 = cntr1 + pnts;
        cntr2 = cntr2 + pnts; 
    end

    %Save to output structure
    if qq==1
        if IRGAFlag
            if no_fw
                %If IRGA Exists @ no finewire
                output.StructParam = [time, Ctson2', Ctsonq', Cq2', r_tsonq'];
                output.StructParamHeader = {'Timestamp', ...
                    [num2str(sensorInfo.Tson(qq, end)), 'm CtSon2: K^2 m^-2/3'],...
                    [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm CtSonq: K kg kg^-1 m^-2/3'],...
                    [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm Cq2: kg^2 kg^-2 m^-2/3'],...
                    [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm r_tSonq'];...
                    [],num2str(sensorInfo.Tson(qq, end)),...
                    num2str(sensorInfo.Tson(qq, end)),...
                    num2str(sensorInfo.irgaH2O(cntrIRGA, end)),...
                    num2str(sensorInfo.irgaH2O(cntrIRGA, end))};
                cntrIRGA = cntrIRGA+1;
            else
                %If IRGA Exists and FW exists
                output.StructParam = [time, Ct2', Ctson2', Ctq', Ctsonq', Cq2', r_tq', r_tsonq'];
                output.StructParamHeader = {'Timestamp', ...
                    [num2str(sensorInfo.Tson(qq, end)), 'm Ctfw2: K^2 m^-2/3'],...
                    [num2str(sensorInfo.Tson(qq, end)), 'm CtSon2: K^2 m^-2/3'],...
                    [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm Ctfwq: K kg kg^-1 m^-2/3'],...
                    [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm CtSonq: K kg kg^-1 m^-2/3'],...
                    [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm Cq2: kg^2 kg^-2 m^-2/3'],...
                    [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm r_tfwq'],...
                    [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm r_tSonq'];...
                    [],num2str(sensorInfo.Tson(qq, end)),...
                    num2str(sensorInfo.Tson(qq, end)),...
                    num2str(sensorInfo.irgaH2O(cntrIRGA, end)),...
                    num2str(sensorInfo.irgaH2O(cntrIRGA, end)),...
                    num2str(sensorInfo.irgaH2O(cntrIRGA, end)),...
                    num2str(sensorInfo.irgaH2O(cntrIRGA, end)),...
                    num2str(sensorInfo.irgaH2O(cntrIRGA, end))};
                cntrIRGA = cntrIRGA+1;
            end
        else
            %If only Sonics
            if no_fw
                output.StructParam = [time, Ctson2'];
                output.StructParamHeader = {'Timestamp', ...
                    [num2str(sensorInfo.Tson(qq, end)), 'm CtSon2: K^2 m^-2/3'];...
                    [], num2str(sensorInfo.Tson(qq, end))};
            else
                output.StructParam = [time, Ct2', Ctson2'];
                output.StructParamHeader = {'Timestamp', ...
                    [num2str(sensorInfo.Tson(qq, end)), 'm Ctfw2: K^2 m^-2/3'],...
                    [num2str(sensorInfo.Tson(qq, end)), 'm CtSon2: K^2 m^-2/3'];...
                    [],num2str(sensorInfo.Tson(qq, end)),...
                    num2str(sensorInfo.Tson(qq, end))};
            end
        end
    else
        if IRGAFlag
            if no_fw
                %IRGA Exists
                output.StructParam = [output.StructParam, Ctson2', Ctsonq', Cq2', r_tsonq'];
                output.StructParamHeader = [output.StructParamHeader,...
                    {[num2str(sensorInfo.Tson(qq, end)), 'm CtSon2: K^2 m^-2/3'],...
                    [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm CtSonq: K kg kg^-1 m^-2/3'],...
                    [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm Cq2: kg^2 kg^-2 m^-2/3'],...
                    [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm r_tSonq'];...
                    num2str(sensorInfo.Tson(qq, end)),...
                    num2str(sensorInfo.Tson(qq, end)),...
                    num2str(sensorInfo.irgaH2O(cntrIRGA, end)),...
                    num2str(sensorInfo.irgaH2O(cntrIRGA, end))}];
                cntrIRGA = cntrIRGA+1;
            else
                %IRGA Exists
                output.StructParam = [output.StructParam, Ct2', Ctson2', Ctq', Ctsonq', Cq2', r_tq',r_tsonq'];
                output.StructParamHeader = [output.StructParamHeader,...
                    {[num2str(sensorInfo.Tson(qq, end)), 'm Ctfw2: K^2 m^-2/3'],...
                    [num2str(sensorInfo.Tson(qq, end)), 'm CtSon2: K^2 m^-2/3'],...
                    [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm Ctfwq: K kg kg^-1 m^-2/3'],...
                    [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm CtSonq: K kg kg^-1 m^-2/3'],...
                    [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm Cq2: kg^2 kg^-2 m^-2/3'],...
                    [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm r_tfwq'],...
                    [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm r_tSonq'];...
                    num2str(sensorInfo.Tson(qq, end)),...
                    num2str(sensorInfo.Tson(qq, end)),...
                    num2str(sensorInfo.irgaH2O(cntrIRGA, end)),...
                    num2str(sensorInfo.irgaH2O(cntrIRGA, end)),...
                    num2str(sensorInfo.irgaH2O(cntrIRGA, end)),...
                    num2str(sensorInfo.irgaH2O(cntrIRGA, end)),...
                    num2str(sensorInfo.irgaH2O(cntrIRGA, end))}];
                cntrIRGA = cntrIRGA+1;
            end
        else
            %Only Sonics
            if no_fw
                output.StructParam = [output.StructParam, Ctson2'];
                output.StructParamHeader = [output.StructParamHeader, ...
                    {[num2str(sensorInfo.Tson(qq, end)), 'm CtSon2: K^2 m^-2/3'];...
                    num2str(sensorInfo.Tson(qq, end))}];
            else
                output.StructParam = [output.StructParam, Ct2', Ctson2'];
                output.StructParamHeader = [output.StructParamHeader, ...
                    {[num2str(sensorInfo.Tson(qq, end)), 'm Ctfw2: K^2 m^-2/3'],...
                    [num2str(sensorInfo.Tson(qq, end)), 'm CtSon2: K^2 m^-2/3'];...
                    num2str(sensorInfo.Tson(qq, end)),...
                    num2str(sensorInfo.Tson(qq, end))}];
            end
        end
    end

    %Save structure functions
    if info.saveStructFunc
        if exist('Dt2', 'var')
            output.StructFunc{qq}.Dt2 = Dt2;
            lenStructFunc = size(Dt2, 2);
        end
        
        if exist('Dtson2', 'var')
            output.StructFunc{qq}.Dtson2 = Dtson2;
            lenStructFunc = size(Dtson2, 2);
        end
        
        if exist('Dtq', 'var')
            output.StructFunc{qq}.Dtq = Dtq;
            lenStructFunc = size(Dtq, 2);
        end
        
        if exist('Dtsonq', 'var')
            output.StructFunc{qq}.Dtsonq = Dtsonq;
            lenStructFunc = size(Dtsonq, 2);
        end
        
        if exist('Dq2', 'var')
            output.StructFunc{qq}.Dq2 = Dq2;
            lenStructFunc = size(Dq2, 2);
        end
        
        output.StructFuncHeader = {['Structure function for ', num2str(lenStructFunc), ' seprations'], ['Separation is [1:1:', num2str(lenStructFunc), ']./freq.*U'], 'Appropriate timestamp can be found in StructParam field'};
        
    end
    clearvars D* C* r_*
    
end