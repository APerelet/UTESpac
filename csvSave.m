function csvSave(template,output,info)
%try

%Check if save directory exists
outputDir = 'outputCSV';

if ~exist([info.rootFolder, filesep, info.siteFolder, filesep, outputDir], 'dir')
    fprintf('\n\tOutput directory does not exist: Creating output dicrectory.\n');
    mkdir([info.rootFolder, filesep, info.siteFolder, filesep, outputDir]);
end

% find serial time stamp from sensible heat flux
serialTime = output.H(:,1);

% create date vector
date = datevec(serialTime);
dataHeader(1:6) = {'yr','month','day','hour','min','sec'};

% eliminate 7.4 m EC150 for Sage Fall 2IRGA
if strcmp(info.tableNames{2},'Sage_20HZ')
    secondIRGAcols = find(~cellfun(@isempty,strfind(output.Sage_20HZHeader(1,:),'7.4')))
    for ii = 1:numel(secondIRGAcols)
        output.Sage_20HZHeader{1,secondIRGAcols(ii)} = ['ex',output.Sage_20HZHeader{1,secondIRGAcols(ii)}];
    end
end

% find Pressure
for i = 1:length(output.tableNames)
    Ptemplate = template.P;
    Pcol = find(strncmpi(output.([output.tableNames{i},'Header'])(1,:),Ptemplate,9)==1);
    Ptble = i;
    if Pcol
        break
    end
end
if isempty(Pcol)
    P = nan(size(date,1),1);
    dataHeader{7} = 'Pressure';
else
    for ii=1:length(Pcol)
        P(:, ii) = output.(output.tableNames{Ptble})(:,Pcol(ii));
        dataHeader{end+1} = output.([output.tableNames{i},'Header']){1,Pcol(ii)};
    end
end

% find rho and cp
dataHeader(end+1:end+2) = {'rho','cp'};
if isempty(Pcol)
    rho = nan(size(date,1),1);
    cp = nan(size(date,1),1);
else
    rho = output.H(:,2);
    cp = output.H(:,3);
end

% find latent heat flux

if isfield(output,'LHflux')
    
    % find columns
    LvCol = find(cell2mat(cellfun(@(x) ~isempty(find(contains(x, 'Lv'))), output.LHfluxHeader, 'UniformOutput', false)));
    
    for ii=1:length(LvCol)
        % sensor height
        tmp = strfind(output.LHfluxHeader{LvCol(ii)}, 'm ');
        LHheight = output.LHfluxHeader{LvCol(ii)}(1:tmp);
        
        %Non corrected flux
        dataHeader(end+1) = {[LHheight, ' LHflux']};
        LHflux(:, ii) = output.LHflux(:,LvCol(ii)).*output.LHflux(:,LvCol(ii)+2);
        
        %WPL corrected flux
        dataHeader(end+1) = {[LHheight, 'LHfluxWPL']};
        LHfluxWPL(:, ii) =  output.LHflux(:,LvCol(ii)+4);
    end
else
    LHflux = [];
    LHfluxWPL = [];
end

% find HMP info 
% intitialize RH, T, headers and z with NaN to allow for Concatnation.  The first column will be deleted
RH = nan;
zRH = nan;
RHheader = {nan};
T = nan;
Theader = {nan};
zT = nan;
for i = 1:length(output.tableNames)
    RHtemplate = template.RH;
    Ttemplate = template.T;
    RHcol = find(strncmpi(output.([output.tableNames{i},'Header'])(1,:),RHtemplate,length(RHtemplate)-1)==1);
    Tcol = find(strncmpi(output.([output.tableNames{i},'Header'])(1,:),Ttemplate,length(Ttemplate)-1)==1);
    RHtble = i;
    if RHcol
        % break
        % RH
        temp = output.(output.tableNames{RHtble})(:,RHcol);
        RH(1:size(temp,1),end+1:end+size(temp,2)) = temp;
        % RH headers and zRH
        temp = output.([output.tableNames{i},'Header'])(1,RHcol);
        RHheader(end+1:end+length(temp)) = output.([output.tableNames{i},'Header'])(1,RHcol);
        zRH = [zRH cell2mat(output.([output.tableNames{i},'Header'])(2,RHcol))];
        
        % T
        temp = output.(output.tableNames{RHtble})(:,Tcol);
        T(1:size(temp,1),end+1:end+size(temp,2)) = temp;
        % T headers and zT
        temp = output.([output.tableNames{i},'Header'])(1,Tcol);
        Theader(end+1:end+length(temp)) = output.([output.tableNames{i},'Header'])(1,Tcol);
        zT = [zT cell2mat(output.([output.tableNames{i},'Header'])(2,Tcol))];
    end
end
% delete first columns
RH(:,1) = []; T(:,1) = []; RHheader(:,1) = []; Theader(:,1) = []; zT(1) = []; zRH(1) = [];

% sort RH
[~, order] = sort(zRH);
RH = RH(:,order);
RHheader = RHheader(:,order);

% sort T
[~, order] = sort(zT);
T = T(:,order);
Theader = Theader(:,order);

% ensure that matrix is at least 6 columns
if size(RH,2) < 6
    RH(:,6) = nan;
    RHheader{:,6} = 'RH_DNE';
    T(:,6) = nan;
    Theader{:,6} = 'Temp_DNE';
end
dataHeader(end+1:end+size(RHheader,2)) = RHheader;
dataHeader(end+1:end+size(RHheader,2)) = Theader;

% find Tson
% intitialize RH, T, headers and z with NaN to allow for Concatnation.  The first column will be deleted
Tson = nan;
zTson = nan;
TsonHeader = {nan};
for i = 1:length(output.tableNames)
    TsonTemplate = template.Tson;
    TsonCol = find(strncmpi(output.([output.tableNames{i},'Header'])(1,:),TsonTemplate,8)==1);
    TsonTble = i;
    if TsonCol
        % break
        temp = output.(output.tableNames{TsonTble})(:,TsonCol);
        Tson(1:size(temp,1),end+1:end+size(temp,2)) = temp;
        % headers and z
        temp = output.([output.tableNames{i},'Header'])(1,TsonCol);
        TsonHeader(end+1:end+length(temp)) = output.([output.tableNames{i},'Header'])(1,TsonCol);
        zTson = [zTson cell2mat(output.([output.tableNames{i},'Header'])(2,TsonCol))];
    end
end
% delete first columns
Tson(:,1) = []; TsonHeader(:,1) = []; zTson(1) = [];

% sort Tson
[~, order] = sort(zTson); sonOrder = order;
Tson = Tson(:,order);
TsonHeader = TsonHeader(:,order);

if size(Tson,2) < 6
    Tson(:,6) = nan;
    TsonHeader{:,6} = 'T_Sonic_DNE';
end
dataHeader(end+1:end+size(TsonHeader,2)) = TsonHeader;
% conver to deg C
Tson(Tson>273) = Tson(Tson>273)-273.15;


% find Tfw
% intitialize RH, T, headers and z with NaN to allow for Concatnation.  The first column will be deleted
Tfw = nan;
zTfw = nan;
TfwHeader = {nan};
for i = 1:length(output.tableNames)
    TfwTemplate = template.fw;
    TfwCol = find(strncmpi(output.([output.tableNames{i},'Header'])(1,:),TfwTemplate,3)==1);
    TfwTble = i;
    if TfwCol
        % break
        temp = output.(output.tableNames{TfwTble})(:,TfwCol);
        Tfw(1:size(temp,1),end+1:end+size(temp,2)) = temp;
        % headers and z
        temp = output.([output.tableNames{i},'Header'])(1,TfwCol);
        TfwHeader(end+1:end+length(temp)) = output.([output.tableNames{i},'Header'])(1,TfwCol);
        zTfw = [zTfw cell2mat(output.([output.tableNames{i},'Header'])(2,TfwCol))];
    end
end
% delete first columns
Tfw(:,1) = []; TfwHeader(:,1) = []; zTfw(1) = [];

% sort Tfw
[~, order] = sort(zTfw);
Tfw = Tfw(:,order);
TfwHeader = TfwHeader(:,order);

% if ~isempty(Tfw)
%     Tfw = output.(output.tableNames{i})(:,TfwCol);
%     TfwHeader = output.([output.tableNames{i},'Header'])(1,TfwCol);
%     [~, order] = sort(cell2mat(output.([output.tableNames{i},'Header'])(2,TfwCol)));
%     Tfw = Tfw(:,order);
%     TfwHeader = TfwHeader(:,order);    
% else
%     Tfw = nan(size(date,1),6);
%     TfwHeader(1:6) = {'FW_DNE','FW_DNE','FW_DNE','FW_DNE','FW_DNE','FW_DNE'};
% end
if size(Tfw,2) < 6
    Tfw(:,6) = nan;
    TfwHeader{:,6} = 'FW_DNE';
end
dataHeader(end+1:end+size(Tfw,2)) = TfwHeader;


% find pot temp from fw
for i = 1:length(output.tableNames)
    if isfield(output,'derivedTheader')
        THfwCol = ~cellfun(@isempty,strfind(output.derivedTheader,'theta_fw'));
    else
        THfwCol = [];
    end
    if (THfwCol)
        break
    end
end
if sum(THfwCol)
    THfw = output.derivedT(:,THfwCol);
    THfwHeader = output.derivedTheader(1,THfwCol);
    THfw = THfw(:,sonOrder);
    THfwHeader = THfwHeader(:,sonOrder);
else
    THfw = nan(size(date,1),6);
    THfwHeader(1:6) = {'theta_fw_DNE','theta_fw_DNE','theta_fw_DNE','theta_fw_DNE','theta_fw_DNE','theta_fw_DNE'};
end
if size(THfw,2) < 6
    THfw(:,6) = nan;
    THfwHeader{:,6} = 'theta_fw_DNE';
end
dataHeader(end+1:end+6) = THfwHeader;

% find direction and speed
for i = 1:length(output.tableNames)
    dirCol = ~cellfun(@isempty,strfind(output.spdAndDirHeader,'direction'));
    spdCol = ~cellfun(@isempty,strfind(output.spdAndDirHeader,'speed'));
    if sum(dirCol)
        break
    end
end
if sum(dirCol)
    dir = output.spdAndDir(:,dirCol);
    dirHeader = output.spdAndDirHeader(1,dirCol);
    spd = output.spdAndDir(:,spdCol);
    spdHeader = output.spdAndDirHeader(1,spdCol);
    dir = dir(:,sonOrder);
    dirHeader = dirHeader(:,sonOrder);
    spd = spd(:,sonOrder);
    spdHeader = spdHeader(:,sonOrder);
else
    dir = nan(size(date,1),6);
    dirHeader(1:6) = {'direction_DNE','direction_DNE','direction_DNE','direction_DNE','direction_DNE','direction_DNE'};
    spd = nan(size(date,1),6);
    spdHeader(1:6) = {'speed_DNE','speed_DNE','speed_DNE','speed_DNE','speed_DNE','speed_DNE'};
end
if size(spd,2) < 6
    dir(:,6) = nan;
    dirHeader{:,6} = 'direction_DNE';
    spd(:,6) = nan;
    spdHeader{:,6} = 'speed_DNE';
end
dataHeader(end+1:end+size(dir,2)) = dirHeader;
dataHeader(end+1:end+size(spd,2)) = spdHeader;

% find w
% intitialize RH, T, headers and z with NaN to allow for Concatnation.  The first column will be deleted
w = nan;
zw = nan;
wHeader = {nan};
for i = 1:length(output.tableNames)
    wTemplate = template.w;
    wCol = find(strncmpi(output.([output.tableNames{i},'Header'])(1,:),wTemplate,3)==1);
    wTble = i;
    if wCol
        % break
        temp = output.(output.tableNames{wTble})(:,wCol);
        w(1:size(temp,1),end+1:end+size(temp,2)) = temp;
        % headers and z
        temp = output.([output.tableNames{i},'Header'])(1,wCol);
        wHeader(end+1:end+length(temp)) = output.([output.tableNames{i},'Header'])(1,wCol);
        zw = [zw cell2mat(output.([output.tableNames{i},'Header'])(2,wCol))];
    end
end

% delete first columns
w(:,1) = []; wHeader(:,1) = []; zw(1) = [];

% sort w
[~, order] = sort(zw); 
w = w(:,order);
wHeader = wHeader(:,order);


if size(w,2) < 6
    w(:,6) = nan;
    wHeader{:,6} = 'Uz_DNE';
end
dataHeader(end+1:end+size(w,2)) = wHeader;

% find sensible heat fluxes
for i = 1:length(output.tableNames)
    h1Col = ~cellfun(@isempty,strfind(output.Hheader,'son:Ts''w'''));
    h2Col = ~cellfun(@isempty,strfind(output.Hheader,'son:Ts''wPF'''));
    h3Col = ~cellfun(@isempty,strfind(output.Hheader,'fw:T''wPF'''));
    if Pcol
        h4Col = ~cellfun(@isempty,strfind(output.Hheader,'fw:Th''wPF'''));
    else
        h4Col = [];
    end
    if sum(h1Col)
        break
    end
end
if sum(h1Col)
    h1 = output.H(:,h1Col);
    h1Header = output.Hheader(1,h1Col);
    h1 = h1(:,sonOrder);
    h1Header = h1Header(:,sonOrder);
else
    h1 = nan(size(date,1),6);
    h1Header = {'son:Ts''w''_DNE','son:Ts''w''_DNE','son:Ts''w''_DNE','son:Ts''w''_DNE','son:Ts''w''_DNE','son:Ts''w''_DNE'};
end
if sum(h2Col)
    h2 = output.H(:,h2Col);
    h2Header = output.Hheader(1,h2Col);
    h2 = h2(:,sonOrder);
    h2Header = h2Header(:,sonOrder);
else
    h2 = nan(size(date,1),6);
    h2Header = {'son:Ts''wPF''_DNE','son:Ts''wPF''','son:Ts''wPF''','son:Ts''wPF''','son:Ts''wPF''','son:Ts''wPF'''};
end
if sum(h3Col)
    h3 = output.H(:,h3Col);
    h3Header = output.Hheader(1,h3Col);
    h3 = h3(:,sonOrder);
    h3Header = h3Header(:,sonOrder);
else
    h3 = nan(size(date,1),6);
    h3Header = {'fw:T''wPF''_DNE','fw:T''wPF''_DNE','fw:T''wPF''_DNE','fw:T''wPF''_DNE','fw:T''wPF''_DNE','fw:T''wPF''_DNE'};
end
if sum(h4Col)
    h4 = output.H(:,h4Col);
    h4Header = output.Hheader(1,h4Col);
    h4 = h4(:,sonOrder);
    h4Header = h4Header(:,sonOrder);
else
    h4 = nan(size(date,1),6);
    h4Header = {'fw:Th''wPF''','fw:Th''wPF''','fw:Th''wPF''','fw:Th''wPF''','fw:Th''wPF''','fw:Th''wPF'''};
end
if size(h1,2) < 6
    h1(:,6) = nan;
    h1Header{:,6} = 'son:Ts''w''_DNE';
    h2(:,6) = nan;
    h2Header{:,6} = 'son:Ts''wPF''';
    h3(:,6) = nan;
    h3Header{:,6} = 'fw:T''wPF''_DNE';
    h4(:,6) = nan;
    h4Header{:,6} = 'fw:Th''wPF''';
end
dataHeader(end+1:end+size(h1,2)) = h1Header;
dataHeader(end+1:end+size(h2,2)) = h2Header;
dataHeader(end+1:end+size(h3,2)) = h3Header;
dataHeader(end+1:end+size(h4,2)) = h4Header;

% find tke
for i = 1:length(output.tableNames)
    tkeCol = ~cellfun(@isempty,strfind(output.tkeHeader,'0.5'));
    if sum(tkeCol)
        break
    end
end
if sum(tkeCol)
    tke = output.tke(:,tkeCol);
    tkeHeader = output.tkeHeader(1,tkeCol);
    tke = tke(:,sonOrder);
    tkeHeader = tkeHeader(:,sonOrder);
else
    tke = nan(size(date,1),6);
    tkeHeader(1:6) = {'0.5(u''^2+v''^2+w''^2)_DNE','0.5(u''^2+v''^2+w''^2)_DNE','0.5(u''^2+v''^2+w''^2)_DNE',...
        '0.5(u''^2+v''^2+w''^2)_DNE','0.5(u''^2+v''^2+w''^2)_DNE','0.5(u''^2+v''^2+w''^2)_DNE'};
end
if size(tke,2) < 6
    tke(:,6) = nan;
    tkeHeader{:,6} = '0.5(u''^2+v''^2+w''^2)_DNE';
end
dataHeader(end+1:end+size(tke,2)) = tkeHeader;

% find momentum fluxes
for i = 1:length(output.tableNames)
    tau1Col = ~cellfun(@isempty,strfind(output.tauHeader,'sqrt(u''w'''));
    tau2Col = ~cellfun(@isempty,strfind(output.tauHeader,':uPF''wPF'''));
    if sum(tau1Col)
        break
    end
end
if sum(tau1Col)
    tau1 = output.tau(:,tau1Col);
    tau1Header = output.tauHeader(1,tau1Col);
    tau1 = tau1(:,sonOrder);
    tau1Header = tau1Header(:,sonOrder);
else
    tau1 = nan(size(date,1),6);
    tau1Header(1:6) = {'sqrt(u''w''^2+v''w''^2)_DNE','sqrt(u''w''^2+v''w''^2)_DNE','sqrt(u''w''^2+v''w''^2)_DNE',...
        'sqrt(u''w''^2+v''w''^2)_DNE','sqrt(u''w''^2+v''w''^2)_DNE','sqrt(u''w''^2+v''w''^2)_DNE'};
end
if size(tau1,2) < 6
    tau1(:,6) = nan;
    tau1Header{:,6} = 'sqrt(u''w''^2+v''w''^2)_DNE';
end
if sum(tau2Col)
    tau2 = output.tau(:,tau2Col);
    tau2Header = output.tauHeader(1,tau2Col);
    tau2 = tau2(:,sonOrder);
    tau2Header = tau2Header(:,sonOrder);
else
    tau2 = nan(size(date,1),6);
    tau2Header(1:6) = {'uPF''wPF''_DNE','uPF''wPF''_DNE','uPF''wPF''_DNE','uPF''wPF''_DNE','uPF''wPF''_DNE','uPF''wPF''_DNE'};
end
if size(tau2,2) < 6
    tau2(:,6) = nan;
    tau2Header{:,6} = 'uPF''wPF''_DNE';
end
dataHeader(end+1:end+size(tau1,2)) = tau1Header;
dataHeader(end+1:end+size(tau2,2)) = tau2Header;

% create data array
data = [date, P, rho, cp, LHflux, LHfluxWPL, RH, T, Tson, Tfw, THfw, dir, spd, w, h1, h2, h3, h4, tke, tau1 tau2];

% store header
fileName2 = strcat(info.siteFolder(5:end),'_Header.csv');
cell2csv(strcat(info.rootFolder,filesep,info.siteFolder,filesep,outputDir,filesep,fileName2),dataHeader)

% store data in 24 hour blocks
numRowsInDay = 24*60./info.avgPer;
if info.PF.globalCalculation
    PFtype = '_GPF';
else
    PFtype = '_LPF';
end
if info.detrendingFormat
    detrendType = '_linDetrend';
else
    detrendType = '_constDetrend';
end
if size(data,1) > numRowsInDay
    % Day 1
    data1 = data(1:numRowsInDay,:);
    fileName = strcat(info.siteFolder(5:end),'_',num2str(info.avgPer),'minAvg_',datestr(output.H(1,1), 'yyyy_mm_dd'),PFtype,detrendType,'.csv');
    csvwrite(strcat(info.rootFolder,filesep,info.siteFolder,filesep,outputDir,filesep,fileName),data1)
    
    % Day 2
    data2 = data(numRowsInDay+1:end,:);
    date = datestr(datenum(info.date)+1,'yyyy_mm_dd');
    fileName = strcat(info.siteFolder(5:end),'_',num2str(info.avgPer),'minAvg_',datestr(output.H(numRowsInDay+1,1), 'yyyy_mm_dd'),PFtype,detrendType,'.csv');
    csvwrite(strcat(info.rootFolder,filesep,info.siteFolder,filesep,outputDir,filesep,fileName),data2)
else
    fileName = strcat(info.siteFolder(5:end),'_',num2str(info.avgPer),'minAvg_',datestr(output.H(1,1), 'yyyy_mm_dd'),PFtype,detrendType,'.csv');
    csvwrite(strcat(info.rootFolder,filesep,info.siteFolder,filesep,outputDir,filesep,fileName),data)
end

%catch err
%    warning(strcat(err.message,'@ line',num2str(err.stack.line),' Unable to create .csv file from the data'))
%end