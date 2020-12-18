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

function [D_xy, C_xy, r] = crossStruct(X, Y, freq, z, sep, r_diff, varargin)

if length(X)~=length(Y)
    error('X and Y must be vectors of equal length');
end

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
    
    %Conversion from temporal to spatial
        %Bosveld, F. C.: The KNMI Garderen experiment:
          %micrometeorological observations 1988–1989,
          %KNMI, the Netherlands,57 pp., 1999.
    denom = (1-1/9*sigma_u^2/U^2+1/3*sigma_v^2/U^2+1/3*sigma_w^2/U^2);
else
    error('Invalid number of inputs, please see cross_struct HELP');
end
    
pnts = 100;
if pnts>length(X)
    pnts=length(X);
end
D_xy = zeros(pnts,1);

%Calculate structure function
for ii=1:pnts
    X_i = X(1+ii:end);
    X_k = X(1:(end-ii));
    Y_i = Y(1+ii:end);
    Y_k = Y(1:(end-ii));
    D_xy(ii)= nanmean((X_i-X_k).*(Y_i-Y_k))./denom;
end 

%Calculate Structure Parameter

if flag
    r = [1:1:100]./freq.*U;
else
    r = [1:1:100]./freq;
end

%Find r/z = sep
[~, ind] = min(abs(r./z-sep));

rVec = r(ind-r_diff:ind+r_diff);

rVec = rVec(rVec>0);
for qq=1:length(rVec)
    r_ = qq;
    
    if isnan(r_)
        C_xy(qq) = nan;
    end
    if flag
        C_xy(qq) = D_xy(r_)*((r_/freq)*U)^(-2/3);
    else
        C_xy(qq) = D_xy(r_)*((r_/freq))^(-2/3);
    end
end