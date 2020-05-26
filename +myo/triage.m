%% TRIAGE remove waveform outliers

function [X,opt] = triage(X,opt)
%%

load([opt.savePath 'spikes'],'Spk')
load([opt.savePath 'labels'],'Lab')
load([opt.savePath 'templates'],'W')

spkIdx = double(Spk.merge);
label = double(Lab.merge);

%% Setup

% get wave snippets
waveLen = round(opt.Fs*opt.waveformDuration(2));
waveforms = myo.spktrig(X, spkIdx, waveLen);

%% Remove outliers

uql = setdiff(unique(label),0);
isNoise = false(size(label));
for iUn = 1:length(uql)
    lid = find(label==uql(iUn));
    if length(lid) > opt.triageMinCount
        isCore = myo.rmvoutliers(waveforms(:,:,lid),opt);
        isNoise(lid(~isCore)) = true;
    else
        isNoise(lid) = true;
        fprintf('Dropped unit %i (%i spikes)\n',uql(iUn),length(lid))
    end
end

% remove noise clusters
spkIdx(isNoise) = [];
label(isNoise) = [];
waveforms(:,:,isNoise) = [];

%% Post-processing

% sort labels by template norm
[w,label] = myo.wavetemplate(waveforms,label);

%% Save results

Spk.triage = uint32(spkIdx);
Lab.triage = uint16(label);
W.triage = w;

save([opt.savePath,'spikes'],'Spk')
save([opt.savePath,'labels'],'Lab')
save([opt.savePath,'templates'],'W')

%% Visualize

if opt.visualize
    myo.viz(X,opt,'triage')
end