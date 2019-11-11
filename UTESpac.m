% Utah Turbulence in Environmental Studies Process and Analysis Code (UTESpac)
% Created by: Derek Jensen, derek591@gmail.com
% Updates to latest version by: Alexei Perelet, alexei.pere@gmail.com
% Version 5.2
% Version Date: 10 Oct 2019
%
% Steps for use:
%
% 1. Create .csv files with the Card Convert Program in Campbell Scientific's LoggerNet.  To do this, open the
% �Destination File Options� box.  Select �Array Compatible CSV� under File Format.  If more than 3 days of turbulence
% data exists, select the �Use Time� toggle and change the �Time Settings� Interval to 2 days.  Under �File Naming�,
% check �TimeDate Filenames� and leave �Use Day of Year� unchecked.  If multiple binary files exist for each table,
% check �Append to Last File�.  Finally, open �Array Compatible CSV Options�.  Under �Timestamp Options� all boxes
% should be checked except �Midnight is 2400�.  Do not check �Include Array ID�.  Lastly, �Array Datalogger Format� with
% �Hour/Minutes and Seconds� should be selected.  This will create .csv files with the date stored in the first 4
% columns (YYYY DOY HHMM SS) of the data table and the date string stored in the filename.
% �CSV_1341.Playa_1HZ_2012_10_10_0000.dat� and �CSV_1341.Playa_20HZ_2012_10_10_0000.dat� are example file names. UTESpac
% will extract the year, month and day (hours and minutes can be different) from the filenames to identify tables that
% belong together. UTESpac will create structures for each data table from midnight to midnight.  If data is missing
% from tables, they will be populated with NaNs. If LoggerNet is not used to convert the data, the data tables can be
% either comma or tab delimited, there should be no header at the top of the table, columns 1 � 4 must contain the date
% in the following format: [YYYY DOY HHMM SS], finally datestrings must be contained in the file names in the following
% format: �tableName�_YYYY_MM_DD where �tableName is a unique identifier for the table.
%
% 2.  The data tables should be placed in a folder with the keyword �site� placed before the site name.  e.g.
% �sitePlayaFall�, �siteES5Fall�, �siteSagbrushSpring�.  UTESpac can be stored elsewhere.  Within the site folder, 4
% items must be included:
% 1) A folder labelled �output� where the output data is stored.
% 2) siteInfo.m must be completed for the specific site.  siteInfo.m is included in UTESpac.zip
% 3) A header file for each different table.  The syntax
% is �tableName�_header.dat (e.g. �Playa_1HZ_header.dat�, �Playa_20HZ_header.dat�).  Note that �tableName� must be
% consistent with the .csv tableNames created in step 1.  The header file is a single line .dat, comma delimited file
% containing variable names and heights for all columns within the respective data table.  The header file is 3 columns
% shorter than the .csv data file.  This is because UTESpac immediately calculates the serial date numbers from the date
% vectors (columns 1 � 4) contained in the data tables.  The serial dates are stored in column 1 and columns 2 � 4 are
% deleted, thus becoming consistent with the header file.  The easiest way to create the header file is with Card
% Convert.  Create an ASCII T0A5 file, there is no need to run the whole binary file, simply stop the conversion
% immediately and only a few hundred lines will be created.  Open the file in a text editor and delete all lines outside
% of the variable headers (typically line 5).  The variable names within the header and the sensor templates (defined on
% lines 155-169) must be consistent.  The template is used by UTESpac to identify specific sensors in the header.  The
% rules for creating the template and header variable names are:
%   - The template and variable name are the exact same except the sensor height is omitted from the template.  e.g.
%     template = 'Ux_', header variable name = 'Ux_0.5', 'Ux_10'
%   - The sensor height must be placed at the end of the header variable name
%   - All sensors (with exception of solar and battery) need an associated height in meters
%   - Heights within the header variable name at a given tower height need to exactly match. e.g. 'FW_5','Ux_5','RH_5'
% 4)  If a global planar fit is used, a PFinfo structure, containing global planar fit coefficients, will be stored in the site
% folder.  There is no need to do anything with it.  Note: For the Global Planar Fit, there must be 1 and only 1 set of
% 5 minute, local planar fit data.  That is, the global planar fit will fail if there is '5minAvg_LPF_linDetrend' and
% '5minAvg_LPF_constDetrend' in the output folder.  There must be one or the other (it doesn't matter which!).
%
% 3.  Fill out the information section of the code (lines 56 - 116) and run the code.  A full example study is included in
% UTESpac.zip
%% ----------------------------------------------INFORMATION------------------------------------------------------------

close all; clearvars; clc;

% current UTESpac Version
info.UTESpacVersion = '5.2';

% enter root folder where site* folders are located
%info.rootFolder = '/uufs/chpc.utah.edu/common/home/IPAQS-group1/IPAQS19/Travis_Scratch/Data/Sonic_Array';
info.rootFolder = '/uufs/chpc.utah.edu/common/home/IPAQS-group1/IPAQS19/Data/Sonic_Array';

% folder structure between site folder and CSV files
% Example
%   .../siteXYZ/*.dat       -> info.foldstruct = '';
%   .../siteXYZ/CSV/*.dat   -> info.foldstruct = [filesep, 'CSV'];
info.foldStruct = [filesep, 'CSV'];
%info.foldStruct = '';

% Enter regular expression for file for
% fields of <Year>, <Month>, <Day> are required
% <TableName> must match what is specified in siteinfo.m 
info.FileForm = '(?<serial>\d+)[.](?<TableName>\w*)_(?<Year>\d{4})_(?<Month>\d{2})_(?<Day>\d{2})_(?<Hour>\d{2})(?<Minute>\d{2}).dat';

% enter averaging period in minutes.  Must yield an integer when dividied into 60 (e.g. 1, 2, 5, 10, 20, 30)
info.avgPer = 30;

% save QC'd raw tables (1 = yes, 0 = no)
info.saveRawConditionedData = true;

% save structure parameters for temperature and humidity
info.saveStructParams = true;

% save netCDF file
info.saveNetCDF = false;

% save .csv files
info.saveCSV = false;

% enter detrending format ('constant' or 'linear')
info.detrendingFormat = 'linear';

% select 'local' or 'global' planar fit, 'local' computes coeffiecients from local file only, 'global' computes
% user-defined, multi-sector, multi-datebin coefficients from all site data - the sector and datebins are defined
% graphically when the code is executed - for 'global' calculations, all data must first be run with a 'local' planar
% fit and 5-min averaging
% MUST be in folder output5
info.PF.globalCalculation = 'global';

% recalulate global PF coefficients if 'global' calculation is used
info.PF.recalculateGlobalCoefficients = false;

% select averaging period for global PF calculation, if used - local PF calculation runs with average specified in
% info.avgPer
info.PF.avgPer = 30;

% chose max and min wind magintudes for global PF calculation, if used - local PF calculation uses all velocity data
info.PF.globalCalcMaxWind = 12;
info.PF.globalCalcMinWind = 0.5;

% give reference specific humidity if humidity measurement does not exist
info.qRef = 12; % (g/kg)
% ---

% --- DATA CONDITIONING OPTIONS See: Vickers and Mahrt 97 and envsupport.licor.com/help/EddyPro3/Content/Topics/Despiking_Raw_Stat_Screening.htm
% 1. spike removal and spike test
info.spikeTest.maxRuns = 20;  % max number of time to run through despike algorithm
info.spikeTest.windowSizeFraction = 1;  % fraction that is multiplied by info.avgPer to get width of averaging window.  Must be <= 1
info.spikeTest.maxConsecutiveOutliers = 10;  % maximum number of consecutive spikes that will be removed.  Longer segments are considered physical
info.spikeTest.maxPercent = 2; % max acceptable percentage of spikes per averaging period
% number of std's to define a spike
info.spikeTest.spikeDef.u = 3.5;
info.spikeTest.spikeDef.v = 3.5;
info.spikeTest.spikeDef.w = 5;
info.spikeTest.spikeDef.Tson = 3.5;
info.spikeTest.spikeDef.fw = 3.5;
info.spikeTest.spikeDef.irgaCO2 = 3.5;
info.spikeTest.spikeDef.irgaH2O = 3.5;
info.spikeTest.spikeDef.KH2O = 3.5;
info.spikeTest.spikeDef.cup = 3.5;
info.spikeTest.spikeDef.birdSpd = 3.5;
info.spikeTest.spikeDef.otherInstrument = 5;

% 2. absolute limits.  Hard flag if single instance occurs over averaging period
info.absoluteLimitsTest.u = [-50, 50]; % m/s
info.absoluteLimitsTest.v = [-50, 50];
info.absoluteLimitsTest.w = [-10, 10];
info.absoluteLimitsTest.Tson = [-20, 80];  % deg C
info.absoluteLimitsTest.fw = [-20, 80];
info.absoluteLimitsTest.irgaCO2 = [0, 1500]; % mg/m^3
info.absoluteLimitsTest.irgaH2O = [0, 50];  % g/m^3
info.absoluteLimitsTest.KH2O = [0, 50];
info.absoluteLimitsTest.cup = [0 50];
info.absoluteLimitsTest.birdSpd = [0 50];

% 3. wind direction
% enter +/- envelope around around tower orientation to define 'bad wind direction'
info.windDirectionTest.envelopeSize = 20;

% 4. excessive NaN test
info.nanTest.maxPercent = 52;

% 5. diagnostic tests
info.diagnosticTest.H2OminSignal = 0.7;
info.diagnosticTest.CO2minSignal = 0.7;
%info.diagnosticTest.meanGasDiagnosticLimit = 0.1; %This seems arbitraty
info.diagnosticTest.meanGasDiagnosticLimit = 10;
info.diagnosticTest.meanSonicDiagnosticLimit = 256; %This seems arbitrary
info.diagnosticTest.meanLiGasDiagnosticLimit = 220;  % Full strength is 255, less indicates problems.  See manual. 
% ---

% SPECIFY SENSOR TEMPLATES [expected units]
template.u = 'Ux_*'; % sonic u  --   [m/s]
template.v = 'Uy_*'; % sonic v  --   [m/s]
template.w = 'Uz_*'; % sonic w  --   [m/s]
template.Tson = 'T_Sonic_*'; % sonic T  --   [C or K]
template.sonDiagnostic = 'sonic_diag_*'; % sonic diagnostic  --  [-]
template.fw = 'fw_*'; % sonic finewires to be used for Eddy Covariance  --  [C]
template.RH = 'RH_HMP_*'; % slow response relative humidity for virtual temperature calculation  --  [Fract or %]
template.T = 'T_HMP_*'; % slow response temperature  --  [C]
template.P = 'pressure'; % pressure  --  [kPa or mBar]
template.irgaH2O = 'H2O_*'; % for use with Campbell EC150 and IRGASON.  WPL corrections applied  --  [g/m^3]
template.irgaH2OsigStrength = 'H2OSig_*'; % EC150 Signal Strength  --  [-]
template.irgaCO2 = 'CO2_*'; % for use with Campbell EC150 and IRGASON.  WPL corrections applied  --  [mg/m^3]
template.irgaCO2sigStrength = 'CO2Sig_*'; % EC150 Signal Strength  --  [-]
template.irgaGasDiag = 'gas_diag_*'; % EC150 gas (CO2 and H2O) diagnostic, 0-> Okay  --  [-]
template.LiH2O = 'LiH2O_*'; % for use with Licor 7500.  WPL corrections applied  --  [mmol/mol]
template.LiCO2 = 'LiCO2_*'; % for use with Licor 7500.  WPL corrections applied  --  [mmol/mol]
template.LiGasDiag = 'Li_gas_diag_*'; % Li7500 gas (CO2 and H2O) diagnostic >~230 -> Okay  --  [-]
template.KH2O = 'KH2O_H2O_*'; % for use with KH2Os.  WPL and O2 corrections applied  --  [g/m^3]
template.cup = 'cup_*';  % wind speed from cup anemometers  --  [m/s]
template.birdSpd = 'wbSpd_*';  % wind speed from prop anemometer  --  [m/s]
template.birdDir = 'wbDir_*';  % wind direction from vain or prop anemometer  --  [deg]
%% ------------------------------------------------ANALYSIS-------------------------------------------------------------
% find data files and load headers
[headers, dataFiles, tableNames, info] = findFiles(info);

% find instrument locations
[sensorInfo, info] = findInstruments(headers,template,info);

% find global PF coefficients
if strcmp(info.PF.globalCalculation,'global')
    PFinfo = findGlobalPF(info,template,sensorInfo);
else
    PFinfo = [];
end

if info.saveStructParams
    if ~info.saveRawConditionedData
        fprintf('\nsaveRawConditionedData is required to save Structure Parameters.\nSet saveRawConditionedData = true\n');
        info.saveRawConditionedData = true;
    end
end

numFiles = size(dataFiles, 1);
parfor i = 1:numFiles %For some reason parfor sometimes gives random errors. Run with regular for loop at that point
%for i = 1:numFiles
        % load files
        [data, dataInfo] = loadData(dataFiles(i,:, :),i,numFiles,info,tableNames);
        
        % store serial date in column 1 of data tables and delete columns 2-4
        [data, dataInfo, fileDate] = findSerialDate(data, dataInfo, info);
        
        % condition data
        [data, output] = conditionData(data,info,tableNames,template,headers);
        
        % store raw averages in output structure
        output =  avg(data,info,tableNames,output,headers,sensorInfo);
        
        % store wind statistics
        output = windStats(output,sensorInfo,tableNames,info);

        % apply local or global (info.PFrange option) planar fit and yaw rotations(Wilczak et. al, 2000)
        [rotatedSonicData, output, dataInfo] = sonicRotation(output,data,sensorInfo,info,dataInfo,tableNames,PFinfo);

        % find sensible heat flux, momentum flux, latent heat flux, CO2 flux
        [output, raw] = fluxes(data, rotatedSonicData, info, output, sensorInfo,tableNames);
        
        % find spatial structure function and structure params
        if info.saveStructParams
            [output , raw] = UTES_struct_setup(data, info, output, raw, tableNames, sensorInfo); %need to make work for multiple heights in same table
        end

        % save data in the root folder
        output = saveData(info,output,dataInfo,tableNames, raw, template, fileDate);
end
gong();
