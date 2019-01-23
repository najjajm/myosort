%% NOISECOV noise covariance
% Function details
%
% SYNTAX
%   outputs = functiontemplate(inputs, varargin)
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

function Cn = noisecov(xn, Fs, dur, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'NOISECOV';

% validation functions
% isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'xn', @isnumeric)
addRequired(P, 'Fs', @isscalar)
addRequired(P, 'dur', @isscalar)
% addParameter(P, 'parameterName', default, validationFunction)

% clear workspace (parser object retains the data while staying small)
parse(P, xn, Fs, dur, varargin{:});
clear ans varargin

%%

% ensure each block has an even number of samples
blockLen = round(Fs*dur);
blockLen = blockLen + mod(blockLen,2);

covLen = size(xn,1);
if mod(covLen,2) == 1
    xn(end,:) = [];
    covLen = covLen-1;
end

nChannels = size(xn,2);

Cn = zeros(blockLen*nChannels);
for ii = 1:nChannels
    for jj = 1:nChannels
        if jj >= ii
            
            % auto/cross-correlation
            xc = ifft(abs(fft(xn(:,ii),2*covLen)).*abs(fft(xn(:,jj),2*covLen)));
            xc = xc(1:blockLen)/covLen;
            
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


%% Testing

% N_SAMP = 1e4;
% noiseSample = zeros(size(C,1),N_SAMP);
% patchIdx = datasample(find(noiseLen>blockLen),N_SAMP,'Replace',false);
% frame = -blockLen/2:blockLen/2-1;
% for ii = 1:N_SAMP
%     y = zeros(blockLen,nChannels);
%     for jj = 1:nChannels
%         y(:,jj) = v(round(noiseLen(patchIdx(ii))/2) + frame,jj);
%     end
%     noiseSample(:,ii) = reshape(y,blockLen*nChannels,1);
% end
% 
% %%
% 
% U = chol(C);
% 
% %%
% 
% x = 0:3e4;
% z = U*noiseSample;
% zs = sum(z.^2,1);
% zp = histc(zs,x);
% zp = zp/sum(zp);
% y = chi2pdf(x,size(z,1));
% %%
% figure
% plot(x,cumsum(zp),'k');
% hold on
% plot(x,cumsum(y),'r')