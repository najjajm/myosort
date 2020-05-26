%% ALIGNWAVEFORMS align multi-channel waveforms
% Function details
%
% SYNTAX
%   outputs = functiontemplate(inputs, varargin)
%
% REQUIRED INPUTS
%   reqIn (class): description
%
% OPTIONAL INPUTS
%   optIn (class): description
%
% PARAMETER INPUTS
%   'parameterName', <argument class>: description (default: )
%
% OUTPUTS
%   out1 (class): description
%
% EXAMPLE(S) 
%
%
% IMPLEMENTATION
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% SEE ALSO:

% Authors: Najja Marshall
% Emails: njm2149@columbia.edu
% Dated:

function [waveformsAligned, spkIdxAligned, dS] = alignwaveforms(X, Fs, waveforms, spkIdx, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'ALIGNWAVEFORMS';

% validation functions
% isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'X', @isnumeric)
addRequired(P, 'Fs', @isscalar)
addRequired(P, 'waveforms', @isnumeric)
addRequired(P, 'spkIdx', @isnumeric)
% addOptional(P, 'optIn', default, validationFunction)
addParameter(P, 'batchSize', 100, @isnumeric)

% clear workspace (parser object retains the data while staying small)
parse(P, X, Fs, waveforms, spkIdx, varargin{:});
clear ans varargin

%% Create sort hierarchy

% waveform energy
wEner = squeeze(sqrt(sum(waveforms.^2,2)))';

[~,waveLen,nObs] = size(waveforms);
grp = ones(nObs,1);
nBatch = nObs;

while any(nBatch > P.Results.batchSize)
    
    % append new layer
    grp = [grp, zeros(nObs,1)];
    
    % split each previous layer
    uql = setdiff(unique(grp(:,end-1)),0);
    
    for ii = 1:length(uql)
        
        lid = (grp(:,end-1)==uql(ii));
    
        if nnz(lid) > P.Results.batchSize
%             pTest = 1-2*min(nnz(lid)-P.Results.batchSize,P.Results.batchSize)/nnz(lid);
            pTest = 1-2*min((nnz(lid)/2)-1,P.Results.batchSize)/nnz(lid);
            lTmp = myosort.split1D(wEner(lid,:),'maxSplit',1,'pThresh',Inf,'pTest',pTest);
            grp(lid,end) = lTmp + max(grp(:,end));
%         else
%             label(lid,end) = ones(nnz(lid),1) + max(label(:,end));
        end
    end
   
    % update group count
    nBatch = cellfun(@(l) nnz(grp(:,end)==l), num2cell(setdiff(unique(grp(:,end)),0)));    
end

%% Align

layer = size(grp,2);
label = (1:nObs)';
while layer > 0
    
    if layer == 1
        maxDepth = Inf;
    else
        maxDepth = 2;
    end
    
    % group indices
    unqGrp = setdiff(unique(grp(:,layer)),0);
    nGrp = length(unqGrp);
    
    lNew = zeros(size(label));    
    for iGrp = 1:nGrp
        
        grpIdx = grp(:,layer)==unqGrp(iGrp);
        
        % labels in group
        lGrp = label(grpIdx);
        unqGrpLab = unique(lGrp,'stable');
        
        % batch size
        batchSize = length(unqGrpLab);
        
        % waveform templates in group
        templates = cellfun(@(l) mean(waveforms(:,:,label==l),3), num2cell(unqGrpLab), 'uni', false);
        templates = cell2mat(reshape(templates,1,1,length(templates)));
        
        % template cross correlation 
        [wXC,enerRat,optLag,~,lags] = myosort.wavexcorr(permute(templates,[2 1 3]),'signed',false);
        optShift = sign(optLag).*lags(abs(optLag));
        wSim = wXC.*enerRat;
        wSim(1:(1+batchSize):batchSize^2) = -Inf;
        
        % shift sequence
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
        
        % optimal time shift
        dS = zeros(batchSize,1);
        for ii = 2:size(shiftSeq,2)
            dS = dS + optShift(sub2ind(size(optShift),shiftSeq(:,ii),shiftSeq(:,ii-1)));
        end
        dS = dS';
        
        % shift waveforms and spikes in batch
        for ii = 1:batchSize
            lid = (label==unqGrpLab(ii));
            spkIdx(lid) = spkIdx(lid) - dS(ii);
            waveforms(:,:,lid) = myosort.spktrig(X, spkIdx(lid), waveLen);
        end
        
        % group
        lTmp = zeros(batchSize,1);
        uqID = unique(shiftSeq(:,end));
        for ii = 1:length(uqID)
            lTmp(shiftSeq(:,end)==uqID(ii)) = ii;
        end
        
        % update labels
        lNewMax = max(lNew);
        for ii = 1:batchSize
            lid = ismember(label,unqGrpLab(ii));
            lNew(lid) = lTmp(ii) + lNewMax;
        end
    end
    
    % copy ungrouped labels
    lNewMax = max(lNew);
    uqlRem = unique(label(grp(:,layer)==0));
    for ii = 1:length(uqlRem)
        lid = ismember(label,uqlRem(ii));
        lNew(lid) = ii + lNewMax;
    end
    
    % overwrite labels
    label = lNew;
    layer = layer-1;
end

% remove duplicates
spkIdx = unique(spkIdx);
fprintf('Removed %i spikes due to refractory violations (%i remaining)\n',...
    length(Spk.detect)-length(spkIdx), length(spkIdx))

% extend waveforms
waveLen = round(opt.Fs*opt.clusterWaveformDuration);
waveforms = myosort.spktrig(X, spkIdx, waveLen);

% center
[~,com] = max(mean(abs(waveforms./noiseStd),3),[],2);
spkIdx = spkIdx - (1+round(waveLen/2)-round(mean(com)));

%%

waveLen = size(waveforms,2);
nObs = length(spkIdx);

[wXC,enerRat,optLag,~,lags] = myosort.wavexcorr(permute(waveforms,[2 1 3]),'signed',false);
optShift = sign(optLag).*lags(abs(optLag));
wSim = wXC.*enerRat;

% construct shift sequence
shiftSeq = (1:nObs)';
while length(unique(shiftSeq(:,end))) > 1
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
dS = zeros(nObs,1);
for ii = 2:size(shiftSeq,2)
    dS = dS + optShift(sub2ind(size(optShift),shiftSeq(:,ii),shiftSeq(:,ii-1)));
end
dS = dS';

% shift spikes
spkIdxAligned = spkIdx - dS;

% remove duplicates
isDuplicate = [false, diff(spkIdxAligned)==0];
spkIdxAligned(isDuplicate) = [];
wXC = wXC(~isDuplicate,~isDuplicate);
enerRat = enerRat(~isDuplicate,~isDuplicate);

% shift waveforms
waveformsAligned = myosort.spktrig(X, spkIdxAligned, waveLen);

% center
[~,com] = max(mean(mean(myosort.smooth1D(abs(waveformsAligned),Fs,'gau','sd',1e-3,'dim',2),1),3));
spkIdxAligned = spkIdxAligned - (1+waveLen/2-com);
waveformsAligned = myosort.spktrig(X, spkIdxAligned, waveLen);