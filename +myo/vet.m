%% VET test remaining templates

function [X,opt] = vet(X,opt)
%%

load([opt.savePath 'batch_indices'],'batchIdx','isNoiseBatch')
load([opt.savePath 'noise_cov'],'Cinv')
load([opt.savePath 'spikes'],'Spk')
load([opt.savePath 'labels'],'Lab')
load([opt.savePath 'templates'],'W')

spkIdx = Spk.triage;
label = double(Lab.triage);
w = W.triage;

%% Test drops

[~,nChan,nUnit] = size(w);
nSamples = size(X,2);

% unique labels
uql = setdiff(unique(label), 0);

% batch length
batchLen = min(ceil(opt.vetBatchDuration*opt.Fs),nSamples-1);
batchLen = batchLen + mod(batchLen,2);
nBatch = opt.vetRounds;

% results metrics
Res = struct('active',{false(1,nUnit)});

% dynamic text strings
unitTxt = sprintf('%%%ii',length(num2str(nUnit)));

% restrict sample indices to non-noise components
batchIdx = batchIdx(~isNoiseBatch);
validTestIdx = cat(2,batchIdx{:});
validTestIdx(validTestIdx < batchLen/2 | validTestIdx > nSamples-batchLen/2+1) = [];

iTest = 1;
runTest = true;
while runTest
    
    fprintf('Test %i: ',iTest)
    
    nUnit = size(w,3);
    
    % reset metrics
    Res.nSpks = zeros(nBatch,nUnit);
    Res.cv = zeros(nBatch,nUnit);
    Res.maxRate = zeros(nBatch,nUnit);
    Res.resEner = zeros(nChan,nUnit,nBatch);
    Res.misaligned = false(nBatch,nUnit);
    
    % draw time limits at random
    batchLim = datasample(validTestIdx,nBatch,'Replace',false)' + [-batchLen/2, batchLen/2-1];
    
    % test fit
    for ii = 1:nBatch
        [spkBatch,~,resEner] = botm(X(:,batchLim(ii,1):batchLim(ii,2))',opt.Fs,w,Cinv,'refDur',1e-3);
        
        % spike count
        Res.nSpks(ii,:) = cellfun(@length,spkBatch);
        
        % coefficient of variation
        isi = cellfun(@(x) diff(x)/opt.Fs,spkBatch,'uni',false);
        isi = cellfun(@(x) x(x<0.2),isi,'uni',false);
        Res.cv(ii,:) = cellfun(@(x) median(std(x)/mean(x)),isi);
        Res.cv(ii,~isfinite(Res.cv(ii,:))) = 0;
        
        % residual energy
        Res.resEner(:,:,ii) = resEner{1};
        
        % misalignment
        Res.misaligned(ii,:) = (all(resEner{1}>1,1) | any(resEner{1}>opt.vetResidualThreshold,1));
    end
    
    % update whether a unit is active (persistent switch)
    Res.active = Res.active | any(Res.nSpks,1);
    
    dropUnit = [];
    
    % check for irregularly firing units
    medCov = cellfun(@(x) median(x(x>0)), mat2cell(Res.cv,nBatch,ones(1,nUnit)));
    if isempty(dropUnit) && any(medCov > opt.vetCVThreshold)
        [~,dropUnit] = max(medCov);
        reason = 'exceeded CoV threshold';
    end
    
    % check for misaligned units
    if isempty(dropUnit) && any(mean(Res.misaligned,1) > opt.vetMisalignThreshold)
        [~,dropUnit] = max(mean(mean(Res.resEner,3),1));
        reason = 'misaligned';
    end
    
    % check for inactive units
    if isempty(dropUnit) && any(~Res.active)
        [~,dropUnit] = min(sum(Res.nSpks,1));
        reason = 'inactive';
    end
    
    % remove unit from set
    if ~isempty(dropUnit)
        
        fprintf(['Dropped unit ' unitTxt ' (' unitTxt ' units remaining). Reason: %s\n'],uql(dropUnit),nUnit-1,reason)
        
        label(label==uql(dropUnit)) = 0;
        uql(dropUnit) = [];
        w(:,:,dropUnit) = [];
        Res.active(dropUnit) = [];
    else
        runTest = false;
        fprintf('Done\n')
    end
    if isempty(w)
        runTest = false;
    end
    iTest = iTest+1;
end

%% Post-processing

spkIdx(label==0) = [];
label(label==0) = [];

% update labels
lNew = zeros(size(label));
for ii = 1:length(uql)
    lNew(label==uql(ii)) = ii;
end
label = lNew;

%% Save results

Spk.vet = spkIdx;
Lab.vet = uint16(label);
W.vet = w;

save([opt.savePath 'spikes'],'Spk')
save([opt.savePath 'labels'],'Lab')
save([opt.savePath 'templates'],'W')

%% Visualize

if opt.visualize
    myo.viz(X,opt,'vet')
end