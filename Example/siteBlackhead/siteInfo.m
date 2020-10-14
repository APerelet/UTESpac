% C-FOG Blackhead 10m tower and Soil
% Ignore GPS Tables
% C-FOG 2016 @ Westmount

% enter orientation of sonics.  Sonic order: tables sorted alphabetically followed by columns sorted in ascending order
info.sonicOrientation = 60.*ones(1, 3);

% enter manufacturer of SAT.  1 for Campbell, 0 for RMYoung.  RMYoung v = Campbell u!
info.sonicManufact = [1, 1, 1];

% enter orientation of tower relative to sonic head
info.tower = (240).*ones(1, 3);

% tower elevation
info.siteElevation = 5; % (m)

% enter expected table names.  Missing tables will be filled with NaNs to create consistency 
% when multiple output files are concatenated with getData.m
info.tableNames = {'Blackhead_10m_FastResponse', 'Blackhead_10m_Radiation', 'Blackhead_10m_SlowResponse', 'Blackhead_Soil_SHF'};

% enter table names to ignore.
info.TableIgnore = {'Blackhead_10m_GPS_Info', 'Blackhead_10m_GPS_NMEA', 'Blackhead_Soil_GPS_Info', 'Blackhead_Soil_GPS_NMEA'};

% enter table scan frequencies corresponding to tableNames
info.tableScanFrequency = [20, 1/60, 1, 1/60];  %[Hz]

% enter number of columns in each .csv table.  Note that the number of columns in the output structure will 
% be 3 less than the number in the .csv file.  This is because the 4 column date vector is replaced with a Matlab's 
% single-column serial time.  Also, note that View Pro frequently cuts of column 1 (the year!) of the .csv file. 
info.tableNumberOfColumns = [27+3, 20+3, 8+3, 15+3];
