%% GETNOISE get noise
% Identifies simultaneous patch of pure noise across all channels in time 
% series data. Uses findpulses to locate all pulse events, identifies the
% longest quiescent period in the data across all channels, then returns
% a section of the data centered on that noisy patch.
%
% SYNTAX
%   outputs = getnoise(X, Fs, [dur], varargin)
%
% REQUIRED INPUTS
%   X (double): time series array (samples x channels)
%   Fs (scalar): sample frequency in Hz
%
% OPTIONAL INPUTS
%   dur (scalar): maximum duration of noise segment in seconds. If no 
%       quiescent period this long exists in the data, will return the 
%       longest section. (default: 1)
%
% PARAMETER INPUTS
%   'normalized', <logical>: indicates whether the data has been normalized
%       by the standard deviation of the noise. Used by findpulses.
%       (Default: false)
%
% OUTPUTS
%   xNoise (double): noisy segment extracted from the data
%   tLim (numeric): time limits of noise segment in seconds
%
% EXAMPLE(S) 
%
% IMPLEMENTATION
% Other m-files required: SMOOTH1D, FINDPULSES
% Subfunctions: SMOOTH1D, FINDPULSES
% MAT-files required: none
%
% SEE ALSO: FINDPULSES

% Authors: Najja Marshall
% Emails: njm2149@columbia.edu
% Dated: February 2019

function [xNoise, tLim] = getnoise(X, Fs, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'GETNOISE';

% validation functions
isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'X', @(x) isnumeric(x) && isa(x,'double'))
addRequired(P, 'Fs', @(x) isscalarnum(x,0,Inf))
addOptional(P, 'dur', 1, @(x) isscalarnum(x,0,Inf))
addParameter(P, 'thresh', 4, @(x) isscalarnum(x,0,Inf))
addParameter(P, 'normalized', false, @islogical)

% clear workspace (parser object retains the data while staying small)
parse(P, X, Fs, varargin{:});
clear ans varargin

%%

[~, nChannels] = size(X);

% locate spikes
spkIdx = findpulses(X,Fs,'normalized',P.Results.normalized,'thresh',P.Results.thresh);

% smooth spike train
s = cell2mat(cellfun(@(idx) sparse(idx,1,true,size(X,1),1),spkIdx,'uni',false));
sF = smooth1D(double(full(s)), Fs, 'box','wid',5e-3);

% identify noiy patches
isNoise = ~any(sF>0.5,2);

% find state transitions
dx = diff([false;isNoise]);
noiseOn = find(dx == 1);
noiseOff = find(dx == -1);

% identify longest region
cutoffLen = min(length(noiseOn), length(noiseOff));
noiseEdges = [noiseOn(1:cutoffLen), noiseOff(1:cutoffLen)];
noiseLen = diff(noiseEdges,[],2);

if ~any(noiseLen >= Fs)
    warning('Longest noise patch detected: %.3f ms. Estimate may be unreliable',max(noiseLen)/(1e3*Fs))
end

% center window within longest noise patch
[~,maxIdx] = max(noiseLen);
cent = round(mean(noiseEdges(maxIdx,:)));

% return longest noise segment up to specified duration
len = min(noiseLen(maxIdx), Fs*P.Results.dur);
len = len + mod(len,2);
frame = cent+(-len/2:len/2-1);

tLim = frame([1 end])/Fs;

xNoise = cell2mat(cellfun(@(ii) X(frame,ii),num2cell(1:nChannels),'uni',false));