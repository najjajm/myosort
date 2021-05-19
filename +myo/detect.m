%% DETECT detect spikes
% Detects spike events. First, denoises the data [Donoho & Johnstone, 1994].
% Then, uses a sequence of peak finding and smoothing operations to 
% identify spikes with mono- or multiphasic signatures.

function [X,opt] = detect(X,opt)

[nChan,nSamp] = size(X);

% estimate standard deviation of noise
noiseStd = double(median(abs(X)/0.6745,2));

% batch data
nBatches = ceil(nChan*nSamp/opt.maxArray);
batchLim = round(linspace(1,nSamp,1+nBatches));
spkIdx = cell(1,nBatches);

for iBatch = 1:nBatches
    
    % normalize and rectify
    Y = abs(double(X(:,batchLim(iBatch):batchLim(iBatch+1)))./noiseStd);
    
    % threshold
    Y = Y .* (Y > opt.detectThreshold);
    
    % filter
    Y = smooth1D(mean(Y,1),opt.Fs,'gau','dim',2,'sd',5e-4);
    
    % find pulses
    [~,spkIdx{iBatch}] = findpeaks(Y,'MinPeakHeight',1);
    
    spkIdx{iBatch} = spkIdx{iBatch} + batchLim(iBatch)-1;
end

spkIdx = cat(2,spkIdx{:});

fprintf('%i spikes found\n',length(spkIdx))

%% Save

Spk.detect = uint32(spkIdx);
Lab.detect = NaN;
W.detect = NaN;

save([opt.savePath,'noise_std'],'noiseStd')
save([opt.savePath,'spikes'],'Spk')
save([opt.savePath,'labels'],'Lab')
save([opt.savePath,'templates'],'W')