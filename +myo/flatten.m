%% FLATTEN flatten multi-channel waveforms

function y = flatten(x)
y = squeeze(reshape(permute(x,[2 1 3]),[1 (size(x,1)*size(x,2)) size(x,3)]));