%% NOISECOV noise covariance
% Estimates spatiotemporal noise covariance matrix from a segment of 
% quiescent time series data. Uses the method of [Pouzat et al., 2002],
% which leverages the observation that the noise covariance for
% multichannel data is itself a Toeplitz matrix constructed from N blocks
% of Toeplitz matrices where N is the number of channels in the data. Each
% block is computed using the auto-/cross-correlation functions of
% the noise segment within/across channels. These functions are then
% truncated to a specified duration in order to construct blocks of
% arbitrary size. To ensure that the covariance matrix is well conditioned,
% the truncation duration should be much less than the duration of the
% noise segments.
%
% SYNTAX
%   Cn = noisecov(Xn, Fs, dur, varargin)
%
% REQUIRED INPUTS
%   Xn (numeric): array of noisy time series data (channels x samples)
%   Fs (scalar): sample frequency in Hz
%   dur (scalar): duration in seconds of each block in the covariance
%       matrix. (For spike sorting, set to the intended waveform duration)
%
% OPTIONAL INPUTS: none
%
% PARAMETER INPUTS: none
%
% OUTPUTS
%   Cn (numeric): noise covariance matrix
%
% EXAMPLE(S) 
%
% IMPLEMENTATION
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% SEE ALSO:

% Authors: Najja Marshall
% Emails: njm2149@columbia.edu
% Dated: February 2019

function Cn = noisecov(Xn, Fs, dur, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'NOISECOV';

% add required, optional, and parameter-value pair arguments
addRequired(P, 'Xn', @isnumeric)
addRequired(P, 'Fs', @isscalar)
addRequired(P, 'dur', @isscalar)

% clear workspace (parser object retains the data while staying small)
parse(P, Xn, Fs, dur, varargin{:});
clear ans varargin

%%

[nSamples, nChannels] = size(Xn);

% ensure each block has an even number of samples
blockLen = round(Fs*dur);
blockLen = blockLen + mod(blockLen,2);
if mod(nSamples,2) == 1
    Xn(end,:) = [];
    nSamples = nSamples-1;
end

% estimate covariance
Cn = zeros(blockLen*nChannels);
for ii = 1:nChannels
    for jj = 1:nChannels
        if jj >= ii
            
            % auto/cross-correlation
            xc = ifft(abs(fft(Xn(:,ii),2*nSamples)).*abs(fft(Xn(:,jj),2*nSamples)));
            xc = xc(1:blockLen)/nSamples;
            
            % toeplitz correlation matrix
            M = toeplitz(xc);
            
            % add to covariance matrix
            xFrame = (ii-1)*blockLen+(1:blockLen);
            yFrame = (jj-1)*blockLen+(1:blockLen);
            Cn(xFrame,yFrame) = M;
        end
    end
end

% symmetrize
Cn = triu(Cn) + triu(Cn,1)';