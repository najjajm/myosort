%% BATCH
% Split data into batches.

function [X,opt] = batch(X,opt)
%%

load([opt.savePath 'noise_std'],'noiseStd')
load([opt.savePath 'spikes'],'Spk')

spkIdx = double(Spk.detect);

%%

[~,nSamp] = size(X);

% create buffers around each spike
s = zeros(1,nSamp);
s(spkIdx) = 1;
spkBuff = smooth1D(s,opt.Fs,'box','wid',2*opt.waveformDuration(2));
spkBuff = spkBuff>0.5;

% identify batch limits as transition points between spike buffers
batchLim = [1,1+find(abs(diff(spkBuff))>0.5),1+nSamp]';
nBatch = length(batchLim)-1;

% compute batch energy
batchEner = zeros(nBatch,1);
for ii = 1:nBatch
    idx = batchLim(ii):batchLim(ii+1)-1;
    batchEner(ii) = mean(sqrt(mean((double(X(:,idx))./noiseStd).^2,2)));
end

% sort by energy
[~,srtIdx] = sort(batchEner,'descend');

% batch indices
batchIdx = mat2cell([batchLim(1:end-1), batchLim(2:end)-1],ones(nBatch,1),2);
batchIdx = cellfun(@(x) x(1):x(2),batchIdx,'uni',false)';
batchIdx = batchIdx(srtIdx);

% identify noise batches
isNoiseBatch = false(size(batchIdx));
for ii = 1:nBatch
    isNoiseBatch(ii) = ~any(ismember(spkIdx,batchIdx{ii}));
end

%% Save results

save([opt.savePath,'batch_indices'],'batchIdx','isNoiseBatch')