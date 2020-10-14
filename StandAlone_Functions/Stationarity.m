%Signal stationarity test
%Foken, Wichura 1996
%Inputs
    %x_i -> first high frequency measurement
    %x_j -> second high frequency measurement
    %M -> length of subinterval (points) for nonstationarity test
%NOTE:
    %x_i and x_j must be the same length
    %M must be smaller than the length of
%Outputs
    %flag -> false -> stationarity is violated
    %        true -> not violated
    %value -> percent error between subinterval and full interval

function [flag, value] = Stationarity(x_i, x_j, M)
if size(x_i, 2)~=1
    x_i = x_i';
end
if size(x_j, 2)~=1
    x_j = x_j';
end

N = length(x_i);

%limit for Stationarity flag
limit = 0.3;

%check integrity of inputs
if N/M < 4
    warning(['Subinterval M = ', num2str(M), ' is large. Consider use of a smaller interval for better results']);
elseif N/M > 8
    warning(['Subinterval M = ', num2str(M), ' is small. Consider use of a larger interval for better results']);
end

if length(x_i)~=length(x_j)
    error('Vectors x_i and x_j must be the same size');
end

% % 
%perturbations with constant detrend
xiP = x_i-ones(N, 1).*nanmean(x_i);
xjP = x_j-ones(N, 1).*nanmean(x_j);
% % 

%covariance of x_i and x_j for subintervals of length M
cntr1 = 1;
cntr2 = M;
for ii=1:ceil(M/N)
    xixj_M(ii) = 1/(M-1).*nansum(x_i(cntr1:cntr2).*x_j(cntr1:cntr2))-...
        1/(M*(M-1)).*nansum(x_i(cntr1:cntr2)).*nansum(x_j(cntr1:cntr2));

% %     xixj_M(ii) = nanmean(xiP(cntr1:cntr2).*xjP(cntr1:cntr2));
    
    cntr1 = cntr1+M;
    cntr2 = cntr2+M;
end

%mean of subinterval covarainces
xixj_mean = nanmean(xixj_M);

%covariance of interval N
xixj_N = 1/(N-1).*nansum(x_i.*x_j)-1/(N*(N-1)).*nansum(x_i).*nansum(x_j);
% % xixj_N = nanmean(xiP.*xjP);

%check if signal stationary over interval

value = abs((xixj_N-xixj_mean)/xixj_N);

if value>=limit
    flag = false;
else
    flag = true;
end
    
