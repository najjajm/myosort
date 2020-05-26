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

function [pkXC, enerRat, optLag, xc, lags] = wavexcorr(w, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'WAVEXCORR';

% add required, optional, and parameter-value pair arguments
addRequired(P, 'w', @(x) isnumeric(x) && ndims(x)<=3)
addParameter(P, 'target', [], @(x) isempty(x) || isnumeric(x))
addParameter(P, 'fullXC', false, @islogical)
addParameter(P, 'signed', true, @islogical)

% clear workspace (parser object retains the data while staying small)
parse(P, w, varargin{:});
clear ans varargin

%%

if ~isempty(P.Results.target)
    w = cat(3,P.Results.target,w);
end

% template dimensions
[waveLen, nChan, nUnit] = size(w);
lags = (1-waveLen):(waveLen-1);

% pairwise inner product for extended templates
wEner = squeeze(sum(sum(w.^2,1),2));
wNorm = squeeze(sqrt(sum(w.^2,1)))+eps; % ensure non-zero entries
if nnz(size(wNorm)>1)==1
    wNorm = wNorm';
end

% compute FFT of normal and time-reversed templates
len = 2*waveLen-1;
nPad = 2^nextpow2(len);
W = mat2cell(w, waveLen, nChan, ones(1,1,nUnit));
wFT = cellfun(@(x) fft(x,nPad,1),W,'uni',false);
uFT = cellfun(@(x) fft(flipud(x),nPad,1),W,'uni',false);

if isempty(P.Results.target) % pairwise cross-correlations
    
    % cross-correlation
    if P.Results.fullXC
        xc = zeros(nUnit,nUnit,len);
    else
        xc = [];
    end
    [pkXC,optLag,enerRat] = deal(zeros(nUnit));
    
    for ii = 1:nUnit-1
        jj = (1+ii):nUnit;
        y = sum(ifft(wFT{ii}.*cell2mat(uFT(jj)),[],1),2);
        y = y(1:len,:,:)./reshape(sqrt(wEner(ii)*wEner(jj)),1,1,nUnit-ii);
        if P.Results.fullXC
            xc(ii,jj,:,:) = permute(y,[4 3 1 2]);
        end
        if P.Results.signed
            [~,ol] = max(abs(y),[],1);
        else
            [~,ol] = max(y,[],1);
        end
        optLag(ii,jj) = reshape(ol,1,nUnit-ii);
        pkXC(ii,jj) = y(sub2ind(size(y),optLag(ii,jj),ones(1,nUnit-ii),1:nUnit-ii));
        enerRat(ii,jj) = mean(min(wNorm(:,ii),wNorm(:,jj))./max(wNorm(:,ii),wNorm(:,jj)),1);
    end
    
    % symmetrize
    pkXC = triu(pkXC,1) + triu(pkXC,1)';
    optLag = triu(optLag,1) - triu(optLag)';
    enerRat = triu(enerRat,1) + triu(enerRat,1)';
    
    optLag(1:(1+nUnit):nUnit^2) = waveLen;
    
else % cross-correlation relative to target template
    
    xc = sum(ifft(wFT{1}.*cell2mat(uFT(2:end)),[],1),2);
    xc = xc(1:len,:,:)./reshape(sqrt(wEner(1)*wEner(2:end)),1,1,nUnit-1);
    xc = squeeze(xc)';
    if P.Results.signed
        [~,optLag] = max(abs(xc),[],2);
    else
        [~,optLag] = max(xc,[],2);
    end
    pkXC = xc(sub2ind(size(xc),(1:nUnit-1)',optLag));
    
    enerRat = mean(min(wNorm(:,1),wNorm(:,2:end))./max(wNorm(:,1),wNorm(:,2:end)),1)';
end