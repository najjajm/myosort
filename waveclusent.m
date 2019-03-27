function wEnt = waveclusent(waveforms,nPCs)

wEnt = zeros(size(waveforms));
for ii = 1:size(wEnt,2)
    
    % projections of waveforms into PCs
    ws = cell2mat(waveforms(:,ii));
    pcs = pca(ws-mean(ws,1));
    wProj = cellfun(@(x) (x-mean(ws,1))*pcs(:,1:nPCs),waveforms(:,ii),'uni',false);
    
    % covariance matrix
    Cw = cellfun(@(x) (x-mean(x,1))'*(x-mean(x,1))/(size(x,1)-1),wProj,'uni',false);
    
    % cluster entropy
    wEnt(:,ii) = cellfun(@(C) nPCs/2 + nPCs/2*log(2*pi) + 1/2*log(det(C)), Cw);
end
wEnt(wEnt~=real(wEnt)) = Inf;