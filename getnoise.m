%% GETNOISE get noise
% Isolates noisy patches from Ephys data
%
% SYNTAX
%   outputs = getnoise(v, varargin)
%
% REQUIRED INPUTS
%   reqIn <class>: description
%
% OPTIONAL INPUTS
%   optIn <class>: description
%
% PARAMETER INPUTS
%   'parameterName' <class>: description (default: )
%
% OUTPUTS
%   out1 <class>: description
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
% Dated:

function vn = getnoise(v, Fs, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'FUNCTIONTEMPLATE';

% validation functions
% isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'v', @isnumeric)
addRequired(P, 'Fs', @isscalar)
% addOptional(P, 'optIn', default, validationFunction)
% addParameter(P, 'parameterName', default, validationFunction)

% clear workspace (parser object retains the data while staying small)
parse(P, v, Fs, varargin{:});
clear ans varargin

%%

[~, nChannels] = size(v);

% locate spikes
spkIdx = findpulses(v,Fs);

% smooth spike train
s = cell2mat(cellfun(@(idx) sparse(idx,1,true,size(v,1),1),spkIdx,'uni',false));
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

% center covariance window within longest noise patch
[~,maxIdx] = max(noiseLen);
cent = round(mean(noiseEdges(maxIdx,:)));

% use longest noise data (up to 1 sec) for noise estimate
covLen = min(noiseLen(maxIdx),Fs);
covLen = covLen + mod(covLen,2);
covFrame = cent+(-covLen/2:covLen/2-1);

vn = cell2mat(cellfun(@(ii) v(covFrame,ii),num2cell(1:nChannels),'uni',false));