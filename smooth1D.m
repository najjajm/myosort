%% SMOOTH1D - 1 dimensional smoothing
% Smooths an input signal using the specified method
%
% SYNTAX
%   xS = smooth1D(x, Fs, method, [norm], varargin)
%
% REQUIRED INPUTS
%   x (vector double) - input signal
%   Fs (scalar double) - sample rate in Hz
%   method (char) - smoothing method. Options: 'gau', 'exp', 'box'
%
% OPTIONAL INPUTS
%   norm (logical) - indicator to normalize smoothed output by height of
%       the filter's impulse response (default: false)
%
% VARIABLE INPUTS
%   (...,'dim',<integer>) - dimension to apply filter (default: 1)
%   (...,'sd',<scalar>) - standard deviation of Gaussian in sec (default: 25e-3)
%   (...,'widGau',<scalar>) - width of Gaussian in multiples of standard
%       deviations (default: 4)
%   (...,'widBox',<scalar>) - width of boxcar in sec (default: 100e-3)
%
% OUTPUTS
%   xS (vector double) - smoothed input
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

function y = smooth1D(x, Fs, method, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'SMOOTH1D';

% add required, optional, and parameter-value pair arguments
addRequired(P, 'x', @isnumeric)
addRequired(P, 'Fs', @isscalar)
addRequired(P, 'method', @(s) ischar(s) && ismember(s,{'box','gau'}))
addOptional(P, 'norm', false, @islogical)
addParameter(P, 'dim', 1, @(d) isscalar(d) && d==round(d))
addParameter(P, 'sd', 0.025, @isscalar) % standard dev. of Gaussian (sec)
addParameter(P, 'wid', [], @isscalar)

% clear workspace (parser object retains the data while staying small)
parse(P, x, Fs, method, varargin{:});
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

filtDim = P.Results.dim;

% pad x to the length of the impulse response
padLen = min(floor(length(fx)/2), size(x,filtDim));
xPad = permute(double(x),[filtDim,setdiff(1:ndims(x),filtDim)]);
xPad = cat(1,mean(xPad(1:padLen,:,:),1).*ones(padLen,size(xPad,2),size(xPad,3)), xPad,...
            mean(xPad(end-padLen+1:end,:,:),1).*ones(padLen,size(xPad,2),size(xPad,3)));

% FFT
nPad = 2^nextpow2(size(xPad,1));
y = ifft(fft(xPad,nPad).*fft(fx(:),nPad));

% remove padding
y = y(2*padLen + (1:size(x,filtDim)),:,:);

% reshape to input dimensions
shiftOrd = zeros(1,ndims(x));
for ii = 1:ndims(x)
    shiftOrd(ii) = find(size(y)==size(x,ii));
end
y = permute(y,shiftOrd);

% normalize by magnitude of impulse response
if P.Results.norm
    y = y/max(fx);
end