%2nd order(cross)structure function
%Inputs
    %X - first variable
    %Y - second variable (if Y=X then 2nd order structure function of X
    %freq - sampling frequency [Hz]
    %z - measurement height
    %sep - NOT USED
    %r_diff - MUST be set to 0
    %options
        %'temporal' - calculates temporal structure function
        %'spatial' - calculated spatial structure function - Recommend only this option be used
            %NOTE: Needs sigma^2_u, sigma^2_v, sigma^2_w, and horizontal U
            %e.g. cross_struct(X, Y, 'spatial', s_u, s_v, s_w, U)
    
%outputs
    %D_xy - structure function of X and Y
    %C_xy - structure parameter of X and Y in inertial subrange

function [D_xy, C_xy, r_] = crossStruct(X, Y, freq, z, sep, r_diff, varargin)

if length(X)~=length(Y)
    error('X and Y must be vectors of equal length');
end

pnts = 500;
if pnts>length(X)
    pnts=length(X);
end
D_xy = zeros(pnts,1);

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
        C_xy = nan*ones(1, r_diff+1);
        D_xy = nan*D_xy;
        r_ = nan.*(1:1:pnts);
        return;
    end
    
    %Conversion from temporal to spatial
        %Bosveld, F. C.: The KNMI Garderen experiment:
          %micrometeorological observations 1988�1989,
          %KNMI, the Netherlands,57 pp., 1999.
    denom = (1-1/9*sigma_u^2/U^2+1/3*sigma_v^2/U^2+1/3*sigma_w^2/U^2);
else
    error('Invalid number of inputs, please see cross_struct HELP');
end
   

%Calculate structure function
for ii=1:pnts
    X_i = X(1+ii:end);
    X_k = X(1:(end-ii));
    Y_i = Y(1+ii:end);
    Y_k = Y(1:(end-ii));
    D_xy(ii)= nanmean((X_i-X_k).*(Y_i-Y_k))./denom;
end 

% Get separation distances / separation time
if flag
    r_ = [1:1:pnts]./freq.*U;
else
    r_ = [1:1:pnts]./freq;
end

%for jj = 1:pnts
%    C_xy_(jj) = D_xy(jj)*r_(jj)^(-2/3);
%end
%[~, ind] = max(abs(C_xy_));

%[~, ind] = min(abs(diff(log10(D_xy))./diff(log10(r_))'-2/3));

% Find points where Structure function is within 5% of 2/3 power law
% and where the separation distance is less than the measurement height
% ONLY WORKS WITH 'spatial' FLAG
mask = abs((diff(log10(D_xy))./diff(log10(r_))'-2/3))/(2/3)<0.05;
mask(end+1) = false;
mask(or(r_<0.1, r_>z)) = false; 

% Get structure parameter using the inertial subrange mask above
if sum(mask)==0
    C_xy = nan;
    r = nan;
else
    C_xy = median(D_xy(mask), 'omitnan').*median(r_(mask), 'omitnan')^(-2/3);
end

%if or(ind==1, ind==pnts)
%    C_xy = nan;
%    r = nan;
%else
%    C_xy = mean(C_xy_(ind-1:ind+1));
%    r = r_(ind);
%end
% % % 
% % % %Find r/z = sep
% % % [~, ind] = min(abs(r./z-sep));
% % % 
% % % rVec = r(ind-r_diff:ind+r_diff);
% % % 
% % % rVec = rVec(rVec>0);
% % % for qq=1:length(rVec)
% % %     r_ = qq;
% % %     
% % %     if isnan(r_)
% % %         C_xy(qq) = nan;
% % %     end
% % %     if flag
% % %         C_xy(qq) = D_xy(r_)*((r_/freq)*U)^(-2/3);
% % %     else
% % %         C_xy(qq) = D_xy(r_)*((r_/freq))^(-2/3);
% % %     end
% % % end
