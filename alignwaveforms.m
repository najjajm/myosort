%% ALIGNWAVEFORMS align spike waveforms
% Aligns spike waveforms from time series data. By default, shifts spike
% indices to be centered at the point of maximum squared amplitude for each
% waveform, but alignment can also be performed on the smoothed waveform
% envelopes. Since each shifting operation can introduce new maxima to the
% waveform window, alignment is run for several iterations until the
% percentage of properly aligned converges within some tolerance or until
% the number of iterations hits a hard limit. Can optionally remove any 
% spikes that coincide within a refractory period through the course of
% alignment. Removes any remaining misaligned waveforms and spike indices.
%
% SYNTAX
%   [waves, waveIdx] = alignwaveforms(X, Fs, spkIdx, varargin)
%
% REQUIRED INPUTS
%   X (numeric): time series array (channels x samples)
%   Fs (scalar): sample frequency in Hz
%   spkIdx (numeric): spike indices
%
% OPTIONAL INPUTS: none
%
% PARAMETER INPUTS
%   'method', <char>: if 'amp' (default), aligns to the index 
%       corresponding to the largest squared amplitude for each waveform. 
%       If 'env', aligns to the index corresponding to the peak in the 
%       waveform envelope (post-rectification and filtering)
%
%   'waveDur', <scalar>: waveform duration in seconds (default: 3e-3)
%
%   'refDur', <scalar>: refractory duration in seconds (default: 0)
%
%   'tol', <scalar>: tolerance on proportion of waveforms properly aligned. 
%       Must be bounded between 0 and 1. Alignment stops after tolerance is
%       reached (or number of iterations hits limit). (default: 0.99)
%
%   'maxIter', <integer>: limit on number of passes through the data
%       before alignment stops. (default: 10)
%
% OUTPUTS
%   waves (numeric): array of waveform shapes (observations x wave length)
%   spkIdx (numeric): aligned spike indices
%
% IMPLEMENTATION
% Other m-files required: SPKTRIG
% Subfunctions: SPKTRIG
% MAT-files required: none
%
% SEE ALSO: SPIKE, SPKTRIG

% Authors: Najja Marshall
% Emails: njm2149@columbia.edu
% Dated: July 2017

function [waves, spkIdx] = alignwaveforms(X, Fs, spkIdx, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'ALIGNWAVEFORMS';

% validation functions
isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'X', @isnumeric)
addRequired(P, 'Fs', @(x) isscalarnum(x,0,Inf))
addRequired(P, 'spkIdx', @isnumeric)
addParameter(P, 'method', 'amp', @(x) ischar(x) && ismember(x,{'amp','env'}))
addParameter(P, 'waveDur', 3e-3, @(x) isscalarnum(x,0,50e-3))
addParameter(P, 'refDur', 0, @(x) isscalarnum(x,0,10e-3))
addParameter(P, 'tol', 0.99, @(x) isscalarnum(x,0,1))
addParameter(P, 'maxIter', 10, @(x) isscalarnum(x,0,Inf) && x==round(x))

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
[waves,spkIdx] = spktrig(X,num2cell(spkIdx),waveLen);

pAligned = 0;
iter = 0;
while pAligned < P.Results.tol && iter < P.Results.maxIter
    
    switch P.Results.method
        case 'amp'
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
    [waves,spkIdx] = spktrig(X,num2cell(spkIdx),waveLen);
    
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
spkIdx = cell2mat(spkIdx);