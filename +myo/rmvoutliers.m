%% TRIAGE remove waveform outliers

function isCore = rmvoutliers(waveforms,opt)
%%

nObs = size(waveforms,3);

% project waveforms into feature space
[~,wavePC] = pca(myo.flatten(waveforms));
feat = myo.flatten(waveforms)'*wavePC(:,1:opt.clusterFeatures);

% whiten
feat = feat-mean(feat,1);
[featPC,~,eigenvals] = pca(feat);
featWht = feat*featPC*diag(eigenvals)^(-1/2);

% feature bins
BUFF = 1.25;
featLim = [min(featWht,[],1); max(featWht,[],1)];
featBin = zeros(1+opt.triageBins,opt.clusterFeatures);
for ii = 1:opt.clusterFeatures
    featBin(:,ii) = linspace(BUFF*featLim(1,ii),BUFF*featLim(2,ii),1+opt.triageBins);
end

% bin features
bins = cell(1,opt.clusterFeatures);
n = zeros(repmat(opt.triageBins,1,opt.clusterFeatures));
for iDim = 1:opt.clusterFeatures
    [~,~,bins{iDim}] = histcounts(featWht(:,iDim),featBin(:,iDim));
end
binIdx = sub2ind(size(n),bins{:});
unqBinIdx = unique(binIdx);
binCount = cellfun(@(bb) nnz(binIdx==bb),num2cell(unqBinIdx));
n(unqBinIdx) = binCount;

% Gaussian blur over bin counts (generalize to N-D?)
if opt.clusterFeatures == 2
    nFlt = imgaussfilt(n,3);
elseif opt.clusterFeatures == 3
    nFlt = imgaussfilt3(n,3);
end

% approximate observation density
p = nFlt(binIdx);

% remove outliers
[~,srtIdx] = sort(p);
nCut = floor(nObs*opt.triageCutProportion);
isCore = true(nObs,1);
isCore(srtIdx(1:nCut)) = false;