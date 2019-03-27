%% WAVEXCORR waveform cross correlation
% Computes the cross correlation for a set of action potential waveform 
% templates. For a set of multi-channel templates, the cross correlation is
% summed across channels. 
%
% SYNTAX
%   [pkXC, xc, lags] = wavexcorr(w)
%
% REQUIRED INPUTS
%   w (numeric): template of dimensions [wavelength x channels x units]
%
% OPTIONAL INPUTS: none
%
% PARAMETER INPUTS: none
%
% OUTPUTS
%   out1 (class): description
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

function [pkXC, enerRat, xc, lags] = wavexcorr(w, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'WAVEXCORR';

% add required, optional, and parameter-value pair arguments
addRequired(P, 'w', @(x) isnumeric(x) && ndims(x)<=3)

% clear workspace (parser object retains the data while staying small)
parse(P, w);
clear ans varargin

%%

% template dimensions
[waveLen, nChan, nUnit] = size(w);

% pairwise inner product for extended templates
wPad = cat(1,zeros(waveLen-1,nChan,nUnit), w, zeros(waveLen-1,nChan,nUnit));
wPad = reshape(wPad, (2*(waveLen-1)+waveLen)*nChan, nUnit);
wEner = diag(wPad'*wPad);
wProd = sqrt(repmat(wEner,1,nUnit) .* repmat(wEner',nUnit,1));

% compute FFT of normal and time-reversed templates
len = 2*waveLen-1;
nPad = 2^nextpow2(len);
W = mat2cell(w, waveLen, nChan, ones(1,1,nUnit));
wFT = cellfun(@(x) fft(x,nPad,1),W,'uni',false);
uFT = cellfun(@(x) fft(flipud(x),nPad,1),W,'uni',false);

% cross-correlation
xc = zeros(nUnit,nUnit,len);
for ii = 1:nUnit
    for jj = 1:nUnit
        if ii~=jj
            y = sum(ifft(wFT{ii}.*uFT{jj},[],1),2);
            xc(ii,jj,:) = reshape(y(1:len),1,1,len)/wProd(ii,jj);
        end
    end
end

% peak cross-correlation
pkXC = max(xc,[],3);

% ratio of waveform energies
M = cellfun(@(x,y) [x,y], repmat(num2cell(wEner),1,nUnit), repmat(num2cell(wEner)',nUnit,1), 'uni', false);
enerRat = cellfun(@(x) min(x)/max(x), M);

% sample lags
lags = (1-waveLen):(waveLen-1);