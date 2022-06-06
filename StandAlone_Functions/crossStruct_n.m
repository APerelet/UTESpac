%2nd order(cross)structure function
%Inputs
    %X - first variable
    %Y - second variable (if Y=X then 2nd order structure function of X
    %freq - sampling frequency [Hz]
    %options
        %'temporal' - calculates temporal structure function
        %'spatial' - calculated spatial structure function
            %NOTE: Needs sigma^2_u, sigma^2_v, sigma^2_w, and horizontal U
            %e.g. cross_struct(X, Y, 'spatial', s_u, s_v, s_w, U)
    
%outputs
    %D_xy - structure function of X and Y
    %C_xy - structure parameter of X and Y in inertial subrange

function [D_xxy, e, r_] = crossStruct_n(X, Y, freq, z, sep, r_diff, varargin)

if length(X)~=length(Y)
    error('X and Y must be vectors of equal length');
end

pnts = 500;
if pnts>length(X)
    pnts=length(X);
end
D_xxy = zeros(pnts,1);

if nargin==6 || strcmp(varargin{1}, 'temporal')
    %calculates temporal structure function
    denom = 1;
    flag=0;
elseif nargin==11 && strcmp(varargin{1}, 'spatial')
    flag=1;
    sigma_u = varargin{2};
    sigma_v = varargin{3};
    sigma_w = varargin{4};
    U = varargin{5};
    
    % Check if nans
    nanCheck = isnan(sigma_u+sigma_v+sigma_w+U);
    if nanCheck
        e = nan*ones(1, r_diff+1);
        D_xxy = nan*D_xxy;
        r_ = nan.*(1:1:pnts);
        return;
    end
    
    %Conversion from temporal to spatial
        %Bosveld, F. C.: The KNMI Garderen experiment:
          %micrometeorological observations 1988–1989,
          %KNMI, the Netherlands,57 pp., 1999.
    denom = (1-1/9*sigma_u^2/U^2+1/3*sigma_v^2/U^2+1/3*sigma_w^2/U^2);
else
    error('Invalid number of inputs, please see cross_struct HELP');
end
   

%Calculate structure function
for ii=1:pnts
    X1_i = X(1+ii:end);
    X1_k = X(1:(end-ii));
    X2_i = X(1+ii:end);
    X2_k = X(1:(end-ii));
    Y_i = Y(1+ii:end);
    Y_k = Y(1:(end-ii));
    D_xxy(ii)= nanmean((X1_i-X1_k).*(X2_i-X2_k).*(Y_i-Y_k))./denom;
end 

%Calculate Structure Parameter

if flag
    r_ = [1:1:pnts]./freq.*U;
else
    r_ = [1:1:pnts]./freq;
end

% Find points where Structure function is within 5% of 1 power law
mask = abs((diff(log10(D_xxy))./diff(log10(r_))'-1))/(1)<0.05;
mask(end+1) = false;
mask(or(r_<0.1, r_>z)) = false; 

if sum(mask)==0
    e = nan;
    r = nan;
else
    e = -4/5.*median(D_xxy(mask), 'omitnan').*median(r_(mask), 'omitnan')^(-1);
end
