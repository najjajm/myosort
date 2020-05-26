%% TEMPLATEMATCHING estimate final spike times
% Estimates spikes using the Bayes Optimal Template Matching algorithm
% [Franke et al., 2015]

function [X,opt] = templatematching(X,opt)
%% Load data

load([opt.savePath 'spikes'],'Spk')
load([opt.savePath 'noise_cov'],'Cinv')
load([opt.savePath 'labels'],'Lab')
load([opt.savePath 'templates'],'W')

fn = fieldnames(W);
w = W.(fn{end});

%%

[waveLen,nChan,nUnit] = size(w);
BUFF_LEN = 10*opt.Fs;
buffFrame = -BUFF_LEN/2:BUFF_LEN/2-1;

% split data into batches as needed
nBatches = 1;
reBatch = true;
while reBatch
    
    batchLim = round(linspace(1,size(X,2),1+nBatches));
    
    % adjust break points to not overlap with spikes
    for ii = 2:nBatches
        s = findpulses(double(X(:,buffFrame+batchLim(ii))),opt.Fs,...
            'dim',2,'OutputFormat','logical','thresh',floor(min(max(max(abs(w),[],1),[],2))));
        spkBuff = false(1,size(s,2));
        for jj = 1:nChan
            spkBuff = spkBuff | (smooth1D(double(s(jj,:)),opt.Fs,'box','wid',2*waveLen/opt.Fs)>0.5);
        end
        noiseIdx = find(~spkBuff);
        dBreak = abs(noiseIdx-(1+BUFF_LEN/2));
        [~,srtIdx] = sort(dBreak);
        batchLim(ii) = noiseIdx(srtIdx(1))+batchLim(ii)+buffFrame(1)-1;
    end
    
    if max(diff(batchLim)*nUnit) < opt.maxArray
        reBatch = false;
    else
        nBatches = nBatches+1;
    end
end

% get final spike times
batchTxt = sprintf('%%%ii',length(num2str(nBatches)));
spkIdx = cell(nBatches,nUnit);
for ii = 1:nBatches
    
    if ii == 1
        fprintf(['Working: batch ' batchTxt '/' batchTxt],ii,nBatches);
    else
        fprintf(repmat('\b',1,1+2*length(sprintf(batchTxt,0))));
        fprintf([batchTxt '/' batchTxt],ii,length(batchLim)-1);
    end
    
    spkIdx(ii,:) = botm(X(:,batchLim(ii):batchLim(ii+1))',opt.Fs,w,Cinv,'refDur',4e-3);
    spkIdx(ii,:) = cellfun(@(x) x(:)+batchLim(ii)-1,spkIdx(ii,:),'uni',false);
end
fprintf('\n')
spkIdx = mat2cell(spkIdx,nBatches,ones(1,nUnit));
spkIdx = cellfun(@(x) cat(1,x{:}),spkIdx,'uni',false);
spkIdx = cellfun(@(x) x(:)',spkIdx,'uni',false);

%% Save results

Spk.templatematching = cell2mat(spkIdx);
Lab.templatematching = cell2mat(cellfun(@(s,l) l*ones(size(s)), spkIdx, num2cell(1:length(spkIdx)), 'uni',false));
W.templatematching = NaN;

save([opt.savePath 'labels'],'Lab')
save([opt.savePath,'spikes'],'Spk')