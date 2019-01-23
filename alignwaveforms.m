%% ALIGNWAVEFORMS - align spike waveforms
% Aligns spike waveforms from time series data
%
% SYNTAX
%   [waves, waveIdx] = alignwaveforms(X, Fs, spkIdx, varargin)
%
% INPUTS
%   X (vector double) - time series data
%   Fs (scalar) - sample frequency in Hz
%   spkIdx (vector double) - spike locations
%
% VARIABLE INPUTS
%   (...,'wdur',waveDur) - waveform duration in ms (default: 5)
%
% OUTPUTS
%   spkIdx (vector double) - spike indices
%   waves (cell array) - waveform shapes
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

function [waves, spkIdx] = alignwaveforms(X, Fs, spkIdx, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'GETSPIKEWAVES';

% validation functions
isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'X', @isnumeric);
addRequired(P, 'Fs', @isscalar);
addRequired(P, 'spkIdx', @isnumeric);
addParameter(P, 'method', 'norm', @(x) ischar(x) && ismember(x,{'norm','env'}))
addParameter(P, 'waveDur', 3e-3, @(x) isscalarnum(x,0,50e-3));   % waveform duration (sec)
addParameter(P, 'refDur', 0, @(x) isscalarnum(x,0,10e-3));    % refractory duration (sec)
addParameter(P, 'tol', 0.99, @(x) isscalarnum(x,0,1)); % tolerance (percent waveforms aligned)
addParameter(P, 'maxIter', 10, @(x) isscalarnum(x,1,Inf))

% clear workspace (parser object retains the data while staying small)
parse(P, X, Fs, spkIdx, varargin{:});
clear ans varargin


%% Extract waveforms, centered on spike locations

% orient spikes vertically
spkIdx = spkIdx(:);

% length of refractory period in samples
refLen = round(Fs * P.Results.refDur);

% remove refractory violations
spkIdx([false; diff(spkIdx) < refLen]) = [];

% waveform window
waveLen = round(Fs * P.Results.waveDur);
waveLen = waveLen + mod(waveLen,2);


%% Center waveforms and remove putative overlaps

if size(X,1)>size(X,2)
    X = X';
end

% extract initial waveforms
[waves,spkIdx] = sta(X,num2cell(spkIdx),waveLen);

pAligned = 0;
iter = 0;
while pAligned < P.Results.tol && iter < P.Results.maxIter
    
    % index of max norm
    switch P.Results.method
        case 'norm'
            [~,maxLoc] = cellfun(@(x) max(double(x).^2), waves, 'uni', false);
                
        case 'env'
            waves = smooth1D(abs(cell2mat(waves)),Fs,'gau','sd',5e-4,'dim',2);
            waves = mat2cell(waves,ones(size(waves,1),1),size(waves,2));
            [~,maxLoc] = cellfun(@(x) max(x), waves, 'uni', false);
    end
    
    % compute distance of max norm from center
    centDist = cellfun(@(ml) ml - waveLen/2, maxLoc, 'uni', false);
    
    % shift spike locations and waveform indices
    spkIdx = cellfun(@(loc,cd) loc + cd, spkIdx, centDist);
    
    % remove refractory violations
    refViol = [false; diff(spkIdx) < refLen];
    spkIdx(refViol) = [];
    maxLoc(refViol) = [];
    
    % update waveforms
    [waves,spkIdx] = sta(X,num2cell(spkIdx),waveLen);
    
    % percent aligned
    pAligned = nnz(cell2mat(maxLoc)==waveLen/2)/length(spkIdx);
    
    iter = iter + 1;
end

%% Post processing

% remove remaining misaligned waveforms
[~,maxLoc] = cellfun(@(x) max(double(x).^2), waves);
isMisaligned = maxLoc~=waveLen/2;
spkIdx(isMisaligned) = [];
waves(isMisaligned) = [];

% output as array
waves = cell2mat(waves);