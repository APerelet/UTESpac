function data_noNans = LinearNanReplace(data)

nanCheck = diff(isnan(data));
nanStart = find(nanCheck == 1);
nanEnd = find(nanCheck == -1);

% If Nans only at end
if isempty(nanEnd)
    nanEnd = length(data)-1;
end

% If nans only at start
if isempty(nanStart)
    nanStart = 1;
end

% If time series ends with nans
if nanStart(end)>nanEnd(end)
    nanEnd(end+1) = length(data)-1;
end

% If time series starts with nans
if nanEnd(1)<nanStart(1)
    nanStart = [1; nanStart];
end

for oo = 1:length(nanStart)
    if isnan(data(nanStart(oo)))
        tmp_ = ones(1, nanEnd(oo)-nanStart(oo)+2).*data(nanEnd(oo)+1);
    elseif isnan(data(nanEnd(oo)+1))
        tmp_ = ones(1, nanEnd(oo)-nanStart(oo)+2).*data(nanStart(oo));
    else
        tmp_ = linspace(data(nanStart(oo)), data(nanEnd(oo)+1), nanEnd(oo)-nanStart(oo)+2);
    end
    data(nanStart(oo):nanEnd(oo)+1) = tmp_;
end

data_noNans = data;