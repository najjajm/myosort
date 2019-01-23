function plotwavetemplate(u,cmapName)

if nargin == 1
    cmapName = 'Spectral';
end

[waveLen,nChan,nUnit] = size(u);

cmap = brewermap(nUnit,cmapName);
cmap = max([0 0 0],cmap-.2);

% cmap = flipud(jet);
% cmap = cmap(round(linspace(1,size(cmap,1),nUnit)),:);
% cmap = max([0 0 0],cmap-.2);

xPos = [0,cumsum(repmat(1.05*size(u,1),1,nUnit-1))];

% vertical offset
yLim = [floor(min(u(:))), ceil(max(u(:)))];
% yLim = [-2660 4920];
yPos = round(flipud([0;cumsum(repmat(1.05*diff(yLim),nChan-1,1))]));

% y-tick position and values
yTickPos = round([yPos+yLim(2)/2, yPos, yPos+yLim(1)/2]);
yTick = yTickPos - yPos;

yTickPos = reshape(yTickPos',3*nChan,1);
yTick = reshape(yTick',3*nChan,1);

figure
hold on
for ch = 1:nChan
    for un = 1:nUnit
        plot(xPos(un)+(1:waveLen),u(:,ch,un)+yPos(ch),'color',cmap(un,:));
    end
end
% title('waveform templates')

% x labels
set(gca,'xtick',xPos+waveLen/2)
unitNo = num2cell(1:nUnit);
set(gca,'xticklabel',cellfun(@(x) sprintf('unit %i',x),unitNo,'uni',false))

% y labels
set(gca,'ytick',flipud(yTickPos))
set(gca,'yTickLabel',cellfun(@num2str,flipud(num2cell(yTick)),'uni',false))
ylabel('noise std.')

% for ch = 1:nChan
%     text(1.025*(waveLen+xPos(end)),yPos(ch),['channel ' num2str(ch)],'fontsize',14)
% end

set(gca,'fontsize',14)

axis tight
set(gca,'xlim',[0, xPos(end)+waveLen])
set(gca,'ylim',[yLim(1), yPos(1)+yLim(2)])