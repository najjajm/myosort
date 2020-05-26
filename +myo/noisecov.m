%% ESTNOISECOV estimate noise covariance matrix
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

function [X,opt] = noisecov(X,opt)
%%

load([opt.savePath 'batch_indices'],'batchIdx')
load([opt.savePath 'templates'],'W')

w = W.triage;

%% Estimate noise covariance from smallest energy batches

% extract noise component from batch indices
batchIdx = cat(2,batchIdx{:});
Xnoise = X(:,batchIdx(end-opt.Fs+1:end))';

% set each block length equal to the template wavelength
blockLen = size(w,1);
[nSamples,nChannels] = size(Xnoise);

% estimate covariance
Cn = zeros(blockLen*nChannels);
for ii = 1:nChannels
    for jj = 1:nChannels
        if jj >= ii
            
            % auto/cross-correlation
            xc = ifft(abs(fft(Xnoise(:,ii),2*nSamples)).*abs(fft(Xnoise(:,jj),2*nSamples)));
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

% invert (blend with identity to ensure non-singular)
Cinv = (opt.noiseAlpha*eye(size(Cn,1)) + opt.noiseAlpha*Cn)^(-1);

%% Save results

save([opt.savePath 'noise_cov'],'Xnoise','Cn','Cinv')