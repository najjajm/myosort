%% MERGETEMPLATES
% Merge templates

function [X,opt] = merge(X,opt)
%%

load([opt.savePath 'spikes'],'Spk')
load([opt.savePath 'labels'],'Lab')
load([opt.savePath 'templates'],'W')

spkIdx = double(Spk.align);
label = double(Lab.cluster);
w = W.cluster;

%% Setup

% get wave snippets
waveLen = round(opt.Fs*opt.waveformDuration(2));
waveforms = myo.spktrig(X, spkIdx, waveLen);

% template similarity
[wXC,~,optLag,~,lags] = wavexcorr(w,'signed',true);
optShift = sign(optLag).*lags(abs(optLag));
    
%% Merge

uql = setdiff(unique(label),0);

while any(wXC(:) > opt.mergeThreshold)
    
    % identify most similar templates
    [unit1,unit2] = ind2sub(size(wXC),find(wXC==max(wXC(:)),1));
    target = min(unit1,unit2);
    source = max(unit1,unit2);
    
    if target == source
        disp('')
    end
    
    % aligned source spikes and waveforms
    spkSource = spkIdx(label==source) - optShift(target,source);
    waveSource = myo.spktrig(X, spkSource, waveLen);
    
    % aggregate target and source waveforms
    u = cell(1,2);
    u{1} = myo.flatten(waveforms(:,:,label==target));
    u{2} = myo.flatten(waveSource);
    u = cellfun(@(x) x-mean(x,1), u, 'uni',false);
    
    % test separability
    [~,pcs] = pca(cat(2,u{:}));
    feat = cat(2,u{:})'*pcs(:,1:opt.clusterFeatures);
    lTmp = split1D(feat,'maxSplit',1);
    
    if (length(unique(lTmp))==1)
        
        spkIdx(label==source) = spkSource;
        waveforms(:,:,label==source) = waveSource;
        label(label==source) = target;
        
        % update list of merged units
        nl = cellfun(@(l) nnz(label==l),num2cell(uql));
        merged = nl==0;
        w(:,:,merged) = 0;
        w(:,:,~merged) = myo.wavetemplate(waveforms,label,false);
        
        % update similarity entries
        [wXCnew,enerRatNew] = wavexcorr(w,'target',w(:,:,target));
        wSimNew = wXCnew.*enerRatNew(:);
        wXC(target,~merged & uql~=target) = wSimNew(~merged & uql~=target);
        wXC(:,target) = wXC(target,:)';
        
        wXC(:,source) = -Inf;
        wXC(source,:) = -Inf;
        
        fprintf('Merged unit %i --> %i\n',source,target);
    else
        wXC(target,source) = -Inf;
        wXC(source,target) = -Inf;
    end
end

%% Post-processing

% sort labels by template norm
[w,label] = myo.wavetemplate(waveforms,label);

%% Save results

Spk.merge = uint32(spkIdx);
Lab.merge = uint16(label);
W.merge = w;

save([opt.savePath,'spikes'],'Spk')
save([opt.savePath,'labels'],'Lab')
save([opt.savePath,'templates'],'W')

%% Visualize

if opt.visualize
    myo.viz(X,opt,'merge')
end