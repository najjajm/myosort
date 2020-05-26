%% CLUSTER cluster waveforms

function [X,opt] = cluster(X,opt)
%%

load([opt.savePath 'spikes'],'Spk')
load([opt.savePath 'labels'],'Lab')
load([opt.savePath 'templates'],'W')

spkIdx = double(Spk.align);

%% Setup

% get wave snippets
waveLen = round(opt.Fs*opt.waveformDuration(1));
waveforms = myo.spktrig(X, spkIdx, waveLen);
nObs = size(waveforms,3);

% wave norm
wNorm = squeeze(sqrt(sum(waveforms.^2,2)))';

%% First pass

label = ones(nObs,1);

nChan = size(X,1);
for iCh = 1:nChan
    
    y = myo.flatten(waveforms(iCh,:,:));
    y = y - mean(y,1);
    
    uql = unique(label);
    for iLa = 1:length(uql)
        
        lid = find(label==uql(iLa));
        
        [~,pcs] = pca(y(:,lid));
        yf = [wNorm(lid,iCh), y(:,lid)'*pcs(:,1:opt.clusterFeatures-1)];
        
        lNew = split1D(yf);
        
        if length(unique(lNew)) > 1
            label(lid) = 0;
            emptyLabel = setdiff(1:max(label)+max(lNew),unique(label));
            for jj = 1:length(unique(lNew))
                label(lid(lNew==jj)) = emptyLabel(jj);
            end
        end
    end
end

%% Second pass

y = myo.flatten(waveforms);
y = y - mean(y,1);

uql = unique(label);
for iLa = 1:length(uql)
    
    lid = find(label==uql(iLa));
    
    [~,pcs] = pca(y(:,lid));
    yf = y(:,lid)'*pcs(:,1:opt.clusterFeatures);
    
    lNew = split1D(yf);
    
    if length(unique(lNew)) > 1
        label(lid) = 0;
        emptyLabel = setdiff(1:max(label)+max(lNew),unique(label));
        for jj = 1:length(unique(lNew))
            label(lid(lNew==jj)) = emptyLabel(jj);
        end
    end
end

%% Post processing

% extend waveforms
waveLen = round(opt.Fs*opt.waveformDuration(2));
waveforms = myo.spktrig(X, spkIdx, waveLen);

% sort labels by template norm
[w,label] = myo.wavetemplate(waveforms,label);

fprintf('%i units detected.\n',size(w,3))

%% Save results

Spk.cluster = NaN;
Lab.cluster = uint16(label');
W.cluster = w;

save([opt.savePath,'spikes'],'Spk')
save([opt.savePath,'labels'],'Lab')
save([opt.savePath,'templates'],'W')

%% Visualize

if opt.visualize
    myo.viz(X,opt,'cluster')
end