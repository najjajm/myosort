%% PLOTWAVES plot waveforms as lines
function plotwaves(waveforms)
%%

nChan = size(waveforms,1);

Y_PAD = 0.05;
Y_TICK = [2/3 1/3];

yLim = squeeze(double(max(max(abs(waveforms),[],2),[],3)));
yPos = flipud([0;cumsum((1+Y_PAD)*2*ones(nChan-1,1))]);

yTickPos = [yPos-Y_TICK, yPos, yPos+fliplr(Y_TICK)];
yTick = round((yTickPos - yPos).*yLim(:));

yTickPos = reshape(fliplr(yTickPos)',(1+2*length(Y_TICK))*nChan,1);
yTick = reshape(fliplr(yTick)',(1+2*length(Y_TICK))*nChan,1);

for iCh = 1:nChan
    plot(squeeze(waveforms(iCh,:,:))/yLim(iCh)+yPos(iCh),'k')
end

set(gca,'ytick',flipud(yTickPos))
set(gca,'yTickLabel',cellfun(@num2str,flipud(num2cell(yTick)),'uni',false))
set(gca,'ylim',[-(1+Y_PAD), yPos(1)+1+Y_PAD])