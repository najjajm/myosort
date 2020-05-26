%% ALIGN temporally align waveforms

function [X,opt] = align(X,opt)
%%

load([opt.savePath 'noise_std'],'noiseStd')
load([opt.savePath 'spikes'],'Spk')
load([opt.savePath 'labels'],'Lab')
load([opt.savePath 'templates'],'W')

spkIdx = double(Spk.detect);

%% Setup

% bound by maximum wave length
waveLen = round(opt.Fs*opt.waveformDuration(2));
spkIdx(spkIdx<2*waveLen | spkIdx>(size(X,2)-2*waveLen)) = [];

% get wave snippets
waveLen = round(opt.Fs*opt.waveformDuration(1));
waveforms = myo.spktrig(X, spkIdx, waveLen);

%% Remove refractory violations

% mean wave norm
wNorm = mean(squeeze(sqrt(sum(waveforms.^2,2)))./noiseStd,1);

% remove smaller waveform from offending pair
refViol = [false, diff(spkIdx)<waveLen];
while any(refViol)
    refViolIdx = find(refViol);
    [~,keepIdx] = max([wNorm(refViolIdx-1);wNorm(refViolIdx)],[],1);
    rmvIdx = zeros(size(refViolIdx));
    rmvIdx(keepIdx==1) = refViolIdx(keepIdx==1);
    rmvIdx(keepIdx==2) = refViolIdx(keepIdx==2)-1;
    spkIdx(rmvIdx) = [];
    waveforms(:,:,rmvIdx) = [];
    wNorm(rmvIdx) = [];
    refViol = [false, diff(spkIdx)<waveLen];
end

%% Align

waveforms = myo.spktrig(X, spkIdx, waveLen);
nObs = size(waveforms,3);

% wave norm, normalized within channels
wNorm = squeeze(sqrt(sum(waveforms.^2,2)));
wNorm = wNorm./max(wNorm,[],2);

% sort spikes and waveforms by energy
[~,srtIdx] = sort(sum(wNorm.*linspace(1,1e3,size(X,1))',1));
spkIdx = spkIdx(srtIdx);
waveforms = waveforms(:,:,srtIdx);

% align in batches
label = (1:nObs)';
doAlign = true;
maxDepth = 2;
while doAlign
    
    uql = unique(label);
    nl = length(uql);
    
    if nl < opt.alignBatchSize
        maxDepth = Inf;
    end
    
    % update batches
    batchLim = [1:opt.alignBatchSize:nl, 1+nl];
    nBatch = length(batchLim)-1;
    
    % append new labels
    lNew = zeros(size(label));
    for iBatch = 1:nBatch
        
        % batch indices
        batchIdx = batchLim(iBatch):batchLim(iBatch+1)-1;
        batchSize = length(batchIdx);
        
        % batch waveforms
        waveBatch = cellfun(@(l) mean(waveforms(:,:,label==l),3), num2cell(batchIdx), 'uni', false);
        waveBatch = cell2mat(reshape(waveBatch,1,1,length(waveBatch)));
        
        % optimal shifts 
        [wXC,enerRat,optLag,~,lags] = wavexcorr(permute(waveBatch,[2 1 3]),'signed',false);
        optShift = sign(optLag).*lags(abs(optLag));
        wSim = wXC.*enerRat;
        wSim(1:(1+batchSize):batchSize^2) = -Inf;
        
        % construct shift sequence
        shiftSeq = (1:batchSize)';
        while size(shiftSeq,2)<maxDepth && length(unique(shiftSeq(:,end))) > 1
            shiftFrom = unique(shiftSeq(:,end));
            shiftTo = shiftFrom;
            [~,bestMatch] = max(wSim(shiftTo,shiftTo),[],2);
            for ii = 1:length(shiftTo)
                shiftTo(ii) = shiftTo(bestMatch(ii));
            end
            [~,mapIdx] = ismember(shiftSeq(:,end),shiftFrom);
            shiftSeq = [shiftSeq, shiftTo(mapIdx)];
        end
        
        % compute optimal shift
        dS = zeros(batchSize,1);
        for ii = 2:size(shiftSeq,2)
            dS = dS + optShift(sub2ind(size(optShift),shiftSeq(:,ii),shiftSeq(:,ii-1)));
        end
        dS = dS';
        
        % shift waveforms and spikes in batch
        for ii = 1:batchSize
            lid = (label==batchIdx(ii));
            spkIdx(lid) = spkIdx(lid) - dS(ii);
            waveforms(:,:,lid) = myo.spktrig(X, spkIdx(lid), waveLen);
        end
        
        % group
        lTmp = zeros(batchSize,1);
        uqID = unique(shiftSeq(:,end),'stable');
        for ii = 1:length(uqID)
            lTmp(shiftSeq(:,end)==uqID(ii)) = ii;
        end
        
        % update labels
        lNewMax = max(lNew);
        for ii = 1:batchSize
            lid = ismember(label,batchIdx(ii));
            lNew(lid) = lTmp(ii) + lNewMax;
        end
    end
    label = lNew;
    
    if maxDepth == Inf
        doAlign = false;
    end
end

% extend waveforms
waveLen = round(opt.Fs*opt.waveformDuration(2));
waveforms = myo.spktrig(X, spkIdx, waveLen);

% center
[~,com] = max(mean(abs(waveforms./noiseStd),3),[],2);
spkIdx = spkIdx - (1+round(waveLen/2)-round(mean(com)));

fprintf('Removed %i spikes due to overlapping waveforms (%i remaining).\n',length(Spk.detect)-length(spkIdx),length(spkIdx))

%% Save results

Spk.align = uint32(spkIdx);
Lab.align = NaN;
W.align = NaN;

save([opt.savePath,'spikes'],'Spk')
save([opt.savePath,'labels'],'Lab')
save([opt.savePath,'templates'],'W')

%% Visualize

if opt.visualize
    myo.viz(X,opt,'align')
end