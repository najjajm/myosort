%% WAVETEMPLATE waveform template
% Constructs the multi-channel waveform template from a set of waveforms,
% sorted by their channel-averaged energies

function [templates, clusSrt] = wavetemplate(waveforms, clus, doSort)
%%

if nargin == 2
    doSort = true;
end

[nChannels,waveLen,~] = size(waveforms);

% unique cluster labels
unqClus = setdiff(unique(clus),0);
nClus = length(unqClus);

% construct mean templates
templates = zeros(waveLen,nChannels,nClus);
for ii = 1:nClus
    templates(:,:,ii) = mean(waveforms(:,:,clus==unqClus(ii)),3)';
end

if doSort
    
    % sort by energy
    ener = squeeze(mean(sum(templates.^2,1),2));
    [~,srtIdx] = sort(ener);
    templates = templates(:,:,srtIdx);
    
    % re-order cluster labels
    clusSrt = zeros(size(clus));
    for ii = 1:length(unqClus)
        clusSrt(clus==unqClus(srtIdx(ii))) = ii;
    end
else
    clusSrt = clus;
end