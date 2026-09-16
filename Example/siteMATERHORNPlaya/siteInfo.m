%MATERHORN Spring
%Playa Tower
% site specific scritp to load site information

% enter orientation of sonics.  Sonic order: tables sorted alphabetically followed by columns sorted in ascending order
info.sonicOrientation = [246 246 246 246 246 250]; %PLAYA

% enter manufacturer of SAT.  1 for Campbell, 0 for RMYoung.  RMYoung v = Campbell u!
info.sonicManufact = [1 1 1 1 1 1]; %Playa

% enter orientation of tower relative to sonic head
info.tower = 47.*ones(1, 6);  % Playa

% tower elevation
info.siteElevation = 1296.6; % (m)

% enter expected table names.  Missing tables will be filled with NaNs to create consistency 
% when multiple output files are concatnated with getData.m
info.tableNames = {'Playa_1HZ','Playa_20HZ'};

% enter table names to ignore.
info.TableIgnore = {''};

% enter table scan frequencies corresponding to tableNames
info.tableScanFrequency = [1, 20];  %[Hz]

% enter number of columns in each .csv table.  Note that the number of columns in the output structure will 
% be 3 less than the number in the .csv file.  This is because the 4 column date vector is replaced with a Matlab's 
% single-column serial time.  Also, note that View Pro frequently cuts of column 1 (the year!) of the .csv file. 
info.tableNumberOfColumns = [17 , 47];