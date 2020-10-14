% C-FOG Blackhead 10m tower and Soil
% Ignore GPS Tables
% C-FOG 2016 @ Westmount

% enter orientation of sonics.  Sonic order: tables sorted alphabetically followed by columns sorted in ascending order
info.sonicOrientation = 270.*ones(1, 1);

% enter manufacturer of SAT.  1 for Campbell, 0 for RMYoung.  RMYoung v = Campbell u!
info.sonicManufact = [1];

% enter orientation of tower relative to sonic head
info.tower = (90).*ones(1);

% tower elevation
info.siteElevation = 5; % (m)

% enter expected table names.  Missing tables will be filled with NaNs to create consistency 
% when multiple output files are concatenated with getData.m
info.tableNames = {'EC150_2'};

% enter table names to ignore.
info.TableIgnore = {''};

% enter table scan frequencies corresponding to tableNames
info.tableScanFrequency = [20];  %[Hz]

% enter number of columns in each .csv table.  Note that the number of columns in the output structure will 
% be 3 less than the number in the .csv file.  This is because the 4 column date vector is replaced with a Matlab's 
% single-column serial time.  Also, note that View Pro frequently cuts of column 1 (the year!) of the .csv file. 
info.tableNumberOfColumns = [13+3];
