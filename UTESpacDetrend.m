% Detrend data for UTESPac based on
% linear, constant or wavelet decomposition

function output =  UTESpacDetrend(data, info)

if strcmp(info.detrendingFormat,'wavelet')
    % Using wavelet decomposition
    [coef,l] = wavedec(data, info.WaveletLevels, 'sym8'); 
    output.coef = coef;
    output.l = l;
    coef(1:l(1)) = 0;
    output.detrendData = waverec(coef, l, 'sym8');
elseif strcmp(info.detrendingFormat,'linear')
    % Using linear detrend
    output.detrendData = detrend(data, 1, 'omitnan');

elseif strcmp(info.detrendingFormat,'constant')
    % using a constant detrend value
    output.detrendData = detrend(data, 0, 'omitnan');

else
    error('Incorrect value for info.detrendingFormat. Must be ''linear'', ''constant'', or ''wavelet''. Check UTESpac.m')
end