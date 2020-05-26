%% Spike-triggered waveform extraction
function [waveforms,spkIdx] = spktrig(X,spkIdx,waveLen)
%%

spkIdx = spkIdx(spkIdx>waveLen/2 & spkIdx<size(X,2)-waveLen/2+1);

nSpks = length(spkIdx);
waveforms = zeros(size(X,1),waveLen,nSpks);

% extract data at each spike
frame = -waveLen/2:waveLen/2-1;
si = num2cell(spkIdx);
for ch = 1:size(X,1)
    y = cellfun(@(idx) X(ch,idx+frame), si,'uni',false);
    waveforms(ch,:,:) = cell2mat(permute(y,[1 3 2]));
end