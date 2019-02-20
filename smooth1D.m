%% SMOOTH1D 1-dimensional smoothing
% Filters a signal via the fast fourier transform. Automatically detects
% and filters along the non-singleton dimension if only one exists;
% otherwise, smooths along any specified input dimension.
%
% SYNTAX
%   Y = smooth1D(X, Fs, method, [norm], varargin)
%
% REQUIRED INPUTS
%   X (numeric): time series array of any size
%   Fs (scalar): sample frequency in Hz
%   method (string): smoothing method. Options: 
%       'gau' (Gaussian filter)
%       'box' (boxcar filter)
%
% OPTIONAL INPUTS
%   norm (logical): if true, normalizes output by amplitude of the filter's
%       impulse response. (default: false)
%
% PARAMETER INPUTS
%   'dim', <integer>: dimension to apply filter. (default: 1, if X has more
%       than one non-singleton dimension)
%
%   'sd', <scalar>: standard deviation of Gaussian in seconds 
%       (default: 0.025)
%
%   'wid' <scalar>: filter widths. Width of Gaussian filter in multiples of
%       standard deviation (default: 4) or width of boxcar filter in
%       seconds (default: 0.1)
%
% OUTPUTS
%   Y (numeric): smoothed input
%
% EXAMPLE(S)
%
%
% IMPLEMENTATION
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% SEE ALSO:

% Authors: Najja Marshall
% Emails: njm2149@columbia.edu
% Dated: July 2017

function Y = smooth1D(X, Fs, method, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'SMOOTH1D';

% validation functions
isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'X', @isnumeric)
addRequired(P, 'Fs', @(x) isscalarnum(x,0,Inf))
addRequired(P, 'method', @(s) ischar(s) && ismember(s,{'box','gau'}))
addOptional(P, 'norm', false, @islogical)
addParameter(P, 'dim', 1, @(x) isscalarnum(x,0,ndims(X)+1) && x==round(x))
addParameter(P, 'sd', 0.025, @(x) isscalarnum(x,0,Inf))
addParameter(P, 'wid', [], @(x) isscalarnum(x,0,Inf))

% clear workspace (parser object retains the data while staying small)
parse(P, X, Fs, method, varargin{:});
clear ans varargin

%% Make impulse response

switch method
    
    case 'box' % Boxcar
        
        % convert parameters to samples
        boxWidSec = P.Results.wid;
        if isempty(boxWidSec)
            boxWidSec = 0.1;
        end
        wid = round(Fs * boxWidSec);
        wid = wid + mod(wid,2);
        
        % sample points
        xx = -wid/2:wid/2;
        
        % impulse response
        fx = ones(size(xx));
        fx = [zeros(1,wid/2), fx, zeros(1,wid/2)];
    
    case 'gau' % Gaussian
        
        % convert parameters to samples
        sd = Fs * P.Results.sd;
        
        gauWidSec = P.Results.wid;
        if isempty(gauWidSec)
            gauWidSec = 4;
        end
        wid = round(sd * gauWidSec);
        
        % Gaussian samples
        xx = -wid:wid;
        
        % impulse response
        fx = 1/(sd*sqrt(2*pi)) * exp(-xx.^2/(2*sd^2));     
end


%% Filter using FFT

% set filter dimension
szX = size(X);
if nnz(szX > 1) == 1
    filtDim = find(szX > 1);
else
    filtDim = P.Results.dim;
end

% pad x to the length of the impulse response
padLen = min(floor(length(fx)/2), szX(filtDim));
Xpad = permute(double(X),[filtDim,setdiff(1:ndims(X),filtDim)]);
Xpad = cat(1,mean(Xpad(1:padLen,:,:),1).*ones(padLen,size(Xpad,2),size(Xpad,3)), Xpad,...
            mean(Xpad(end-padLen+1:end,:,:),1).*ones(padLen,size(Xpad,2),size(Xpad,3)));

% FFT
nPad = 2^nextpow2(size(Xpad,1));
Y = ifft(fft(Xpad,nPad).*fft(fx(:),nPad));

% remove padding
Y = Y(2*padLen + (1:szX(filtDim)),:,:);

% reshape to input dimensions
shiftOrd = zeros(1,ndims(X));
for ii = 1:ndims(X)
    shiftOrd(ii) = find(size(Y)==szX(ii));
end
Y = permute(Y,shiftOrd);

% normalize by magnitude of impulse response
if P.Results.norm
    Y = Y/max(fx);
end