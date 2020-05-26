%% BOTM Bayes optimal template matching
% Implementation of the Bayes optimal template matching algorithm for spike
% sorting [Franke et al., 2015]
%
% SYNTAX
%   outputs = functiontemplate(inputs, varargin)
%
% REQUIRED INPUTS
%   reqIn <class>: description
%
% OPTIONAL INPUTS
%   optIn <class>: description
%
% PARAMETER INPUTS
%   'parameterName' <class>: description (default: )
%
% OUTPUTS
%   out1 <class>: description
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

function [spkIdx,d,resEner] = botm(X, Fs, template, Cinv, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'BOTM';

% validation functions
isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'X', @isnumeric)
addRequired(P, 'Fs', @isnumeric)
addRequired(P, 'template', @isnumeric)
addRequired(P, 'Cinv', @isnumeric)
addParameter(P, 'prior', 0.01, @(x) isscalarnum(x,0,1))
addParameter(P, 'gain', 1, @(x) isscalarnum(x,0,Inf))
addParameter(P, 'sic', true, @islogical)
addParameter(P, 'refDur', 0, @(x) isscalarnum(x,-eps,Inf))
addParameter(P, 'plot', {''}, @(x) iscell(x) && all(ismember(x,{'','fit','fitTight','res','resen','rec','dis'})))
addParameter(P, 'cmap', 'Spectral', @(x) ischar(x) || (isnumeric(x) && size(x,2)==3))
addParameter(P, 'bright', -0.2, @(x) isscalarnum(x,-1,1))
addParameter(P, 'verbose', false, @islogical)
addParameter(P, 'priority', 'stability', @(x) ischar(x) && ismember(x,{'stability','speed'}))
addParameter(P, 'figNo', [], @(x) isempty(x) || isnumeric(x))

% clear workspace (parser object retains the data while staying small)
parse(P, X, Fs, template, Cinv, varargin{:});
clear ans varargin

%% Setup

nDataPoints = size(X,1);
[waveLen, nChan, nUnit] = size(template);

% matched filter
w = cellfun(@(n) reshape(template(:,:,n),waveLen*nChan,1),num2cell(1:nUnit),'uni',false);
u = cellfun(@(x) Cinv*x, w, 'uni', false);

w = cellfun(@(x) reshape(x,waveLen,nChan), w, 'uni', false);
u = cellfun(@(x) reshape(x,waveLen,nChan), u, 'uni', false);

% filter cross-correlation
if P.Results.sic
    len = 2*waveLen-1;
    nPad = 2^nextpow2(len);
    wFT = cellfun(@(x) fft(x,nPad,1),w,'uni',false);
    uFT = cellfun(@(x) fft(flipud(x),nPad,1),u,'uni',false);
    fProj = zeros(nUnit,nUnit,len);
    for ii = 1:nUnit
        for jj = 1:nUnit
            if ii~=jj
                y = sum(ifft(wFT{ii}.*uFT{jj},[],1),2);
                fProj(ii,jj,:) = reshape(y(1:len),1,1,len);
            end
        end
    end
end
lags = (1-waveLen):(waveLen-1);

% whitened filter energy
filtEner = zeros(nChan,nUnit);
for ii = 1:nChan
    filtEner(ii,:) = cellfun(@(a,b) P.Results.gain * -1/2*a(:,ii)'*b(:,ii), w, u);
end

% filter constant
filtConst = filtEner + log(P.Results.prior/nUnit);

% decision threshold
thr = log(1-P.Results.prior);

%% Detect spikes using Fisher's Discriminant

if P.Results.verbose
    fprintf('Computing discriminant. ')
    t0 = tic;
end

% Fourier transform data and time-reversed matched filters
nPad = 2^nextpow2(nDataPoints+waveLen/2-1);
d = zeros(nDataPoints,nUnit);
for ch = 1:nChan
    XFT = fft(double(X(:,ch)),nPad,1);
    for un = 1:nUnit
        y = ifft(XFT .* fft(flipud(u{un}(:,ch)),nPad));
        d(:,un) = d(:,un) + y(waveLen/2-1+(1:nDataPoints)) + filtConst(ch,un);
    end
end

if P.Results.verbose
    fprintf('Runtime: %.2f min\n',toc(t0)/60)
    clear XFT
end

if P.Results.verbose
    fprintf('Detecting spikes. ')
    t0 = tic;
end

% detect spikes
spkIdx = cell(1,nUnit);
warning('off','signal:findpeaks:largeMinPeakHeight')
for un = 1:nUnit
    [~,spkIdx{un}] = findpeaks(d(:,un),'MinPeakHeight',thr);
    spkIdx{un} = spkIdx{un}(spkIdx{un}>waveLen & spkIdx{un}<nDataPoints-waveLen);
end

% extract the discriminant at each spike
% dSpk = (thr-eps) * ones(nDataPoints,nUnit);
% for un = 1:nUnit
%     dSpk(spkIdx{un},un) = d(spkIdx{un},un);
% end

for un = 1:nUnit
    s = false(nDataPoints,1);
    s(spkIdx{un}) = true;
    d(~s,un) = (thr-eps);
end

if P.Results.verbose
    fprintf('Runtime: %.2f min\n',toc(t0)/60)
end

%% Remove refractory period violations

if P.Results.refDur > 0
    
    if P.Results.verbose
        fprintf('Removing refractory period violations. ')
        t0 = tic;
    end
    
    for un = 1:nUnit
        
        % convert spike indices to pulse trains
        s = zeros(1,nDataPoints);
        s(spkIdx{un}) = 1;
    
        % detect putative refractory period violations
        y = double(smooth1D(s,Fs,'box','wid',P.Results.refDur/2,'dim',2));
        y = y.*(y>0.5);
        
        dy = [NaN,diff(y)];
        yOn = find(dy(2:end)>0.5 & dy(1:end-1)==0 & y(1:end-1)<0.5);
        yOff = find(dy(1:end-1)<-0.5 & dy(2:end)==0 & y(2:end)<0.5);
        yLim = mat2cell([yOn(:),yOff(:)],ones(length(yOn),1),2);
        
        refLim = cell2mat(yLim);
        
        rmvSpk = false(size(spkIdx{un}));
        
        for ii = 1:size(refLim,1)
            
            binIdx = spkIdx{un}>=refLim(ii,1) & spkIdx{un}<=refLim(ii,2);
            spksInBin = spkIdx{un}(binIdx);
            si = find(binIdx);
            rmvIdx = false(size(si));
            
            dBin = d(spksInBin,un);
            [~,srtIdx] = sort(dBin,'descend');
            spksInBin = spksInBin(srtIdx);
            si = si(srtIdx);
            
            for jj = 1:length(si)-1
                if ~rmvIdx(jj)
                    dt = abs(spksInBin - spksInBin(jj))/Fs;
                    rmvIdx(dt>0 & dt<=P.Results.refDur) = true;
                end                
            end    
            
            rmvSpk(si(rmvIdx)) = true;
        end
        
        d(spkIdx{un}(rmvSpk),un) = thr-eps;
        spkIdx{un}(rmvSpk) = [];
    end
    
    if P.Results.verbose
        fprintf('Runtime: %.2f min\n',toc(t0)/60)
    end
end
    

%% Subtractive interference cancellation

if P.Results.verbose
    fprintf('Running SIC. ')
    t0 = tic;
end

if P.Results.sic && nnz(~cellfun(@isempty,spkIdx)) > 1
    
    % convert spike indices to pulse trains
    [s,sFilt] = deal(zeros(nUnit,nDataPoints));
    for un = 1:nUnit
        s(un,spkIdx{un}) = 1;
        sFilt(un,:) = double(smooth1D(s(un,:),Fs,'box','wid',waveLen/Fs,'dim',2) > 0.5);
    end
    
    % create spike IDs (index and unit number)
    spkIdx = cellfun(@(x) x(:)',spkIdx,'uni',false);
    spkIdx = spkIdx(:)';
    spkUnit = cellfun(@(n,x) repmat(n,1,length(x)),num2cell(1:nUnit),spkIdx,'uni',false);
    spkID = [cell2mat(spkIdx); cell2mat(spkUnit)];
    
    % detect putative spike overlaps
    y = round(sum(sFilt,1));
    dy = [NaN,diff(y)];
    yOn = find(dy(2:end)>0.5 & dy(1:end-1)==0 & y(1:end-1)<0.5);
    yOff = find(dy(1:end-1)<-0.5 & dy(2:end)==0 & y(2:end)<0.5);
    yLim = mat2cell([yOn(:),yOff(:)],ones(length(yOn),1),2);
    overlapBins = cellfun(@(l) max(y(l(1):l(2)))>1.5, yLim);
    overlapLim = cell2mat(yLim(overlapBins));
    
    if ~isempty(overlapLim)
        
        % assign spikes to overlap bins
        bins = reshape([overlapLim(:,1), 1+overlapLim(:,2)]', 1, 2*size(overlapLim,1));
        [~,b] = histc(spkID(1,:),bins);
        unqB = setdiff(unique(b),0);
        
        % resolve overlaps within bins
        for ii = 1:length(unqB)
            
            spksInBin = spkID(:,b==unqB(ii));
            nSpks = size(spksInBin,2);
            
            % rank units by the amplitude of their discriminant at the spike
            dBin = cellfun(@(x) d(x(1),x(2)),mat2cell(spksInBin,2,ones(1,nSpks)));
            [~,srtIdx] = sort(dBin,'Descend');
            spksInBin = spksInBin(:,srtIdx);
            
            for jj = 1:nSpks-1
                t1 = spksInBin(1,jj);
                u1 = spksInBin(2,jj);
                if d(t1,u1) < thr
                    continue
                end
                for kk = (1+jj):nSpks
                    u2 = spksInBin(2,kk);
                    if u1 == u2
                        continue
                    end
                    t2 = spksInBin(1,kk);
                    tau = t1-t2;
                    if abs(tau) > waveLen-1
                        continue
                    end
                    d(t2,u2) = d(t2,u2) - fProj(u2,u1,lags==tau);
                end
            end
        end
        
        % update spike estimate
        spkIdx = cell(1,nUnit);
        for un = 1:nUnit
            [~,spkIdx{un}] = findpeaks(d(:,un),'MinPeakHeight',thr);
            spkIdx{un} = spkIdx{un}(spkIdx{un}>waveLen & spkIdx{un}<nDataPoints-waveLen);
        end
    end
end

if P.Results.verbose
    fprintf('Runtime: %.2f min\n',toc(t0)/60)
end

%% Residual energy

if nargout > 2 || ismember('resen',P.Results.plot)
    if P.Results.verbose
        fprintf('Computing residual energy. ')
        t0 = tic;
    end
    
    X = double(X);
    
    % identify spike regions
    allSpk = findpulses(X,Fs,'dim',1,'OutputFormat','logical');
    spkRegion = smooth1D(double(allSpk),Fs,'box','wid',waveLen/(2*Fs),'dim',1);
    
    resEner = cell(1,2);
    resEner(:) = {zeros(nUnit,nChan)};
    s = zeros(nDataPoints,1);
    for un = 1:nUnit
        s = s*0;
        s(spkIdx{un}) = 1;
        for ch = 1:nChan
            res = X(:,ch) - conv(s,template(:,ch,un),'same');
            resEner{1}(un,ch) = sum(res.^2);
            resEner{2}(un,ch) = sum(res(spkRegion(:,ch)>0.5).^2);
        end
    end
    resEner{1} = (resEner{1}./sum(X.^2,1))';
    for ch = 1:nChan
        resEner{2}(:,ch) = resEner{2}(:,ch)./sum(X(spkRegion(:,ch)>0.5,ch).^2);
    end
    resEner{2} = resEner{2}';
    
    if P.Results.verbose
        fprintf('Runtime: %.2f min\n',toc(t0)/60)
    end
else
    resEner = [];
end

%% Plotting

if all(cellfun(@isempty,P.Results.plot))
    return
end

% color map
cmap = P.Results.cmap;
if ischar(cmap)
    cmap = brewermap(nUnit,cmap);
    cmap = min([1 1 1],max([0 0 0], cmap + P.Results.bright));
end

% axes
figNo = P.Results.figNo;
nPlots = length(P.Results.plot);
if isempty(figNo) || length(figNo)~=nPlots
    holdFig = false;
else
    holdFig = true;
end
figIdx = 1;

t = (1:nDataPoints)/Fs;

% vertical offsets
Y_PAD = 0.1;
Y_TICK = [2/3 1/3];

yLim = double(max(abs(X),[],1));
yPos = flipud([0;cumsum((1+Y_PAD)*2*ones(nChan-1,1))]);

yTickPos = [yPos-Y_TICK, yPos, yPos+fliplr(Y_TICK)];
yTick = round((yTickPos - yPos).*yLim(:));

yTickPos = reshape(fliplr(yTickPos)',(1+2*length(Y_TICK))*nChan,1);
yTick = reshape(fliplr(yTick)',(1+2*length(Y_TICK))*nChan,1);

% discriminant
if ismember('dis',P.Results.plot)
    if holdFig
        figure(figNo(figIdx))
        clf
        figIdx = figIdx+1;
    else
        figure
    end
    ax = [0 0];
    ax(1) = subplot(211);
    hold on
    for ii = 1:nUnit
        plot(d(:,ii),'color',cmap(ii,:));
    end
    set(gca,'xlim',[1 size(d,1)])
    plot(get(gca,'xlim'),thr*[1 1],'color',0.8*[1 1 1])
    box off
    title(sprintf('Fisher''s Discriminant\nraw'))
    legend(cellfun(@(n) sprintf('unit %i',n),num2cell(1:nUnit),'uni',false))
    ax(2) = subplot(212);
    hold on
    for ii = 1:nUnit
        plot(d(:,ii),'color',cmap(ii,:),'linewidth',1.5)
    end
    title('post-SIC')
    set(gca,'xlim',[1 size(d,1)])
    plot(get(gca,'xlim'),thr*[1 1],'color',0.8*[1 1 1])
    box off
    title(sprintf('Fisher''s Discriminant\nthresholded'))
    linkaxes(ax,'x')
end

% single unit fits
if ismember('fit',P.Results.plot) || ismember('fitTight',P.Results.plot)  
    if holdFig
        figure(figNo(figIdx))
        clf
        figIdx = figIdx+1;
    else
        figure
    end
    hold on
    s = zeros(nDataPoints,1);
%     ax = [];
    fh = [];
    for ch = 1:nChan
        
%         ax(ch) = subplot(nChan,1,ch);
        plot(t,double(X(:,ch))/yLim(ch)+yPos(ch),'color',[0 0 0 0.25],'linewidth',2.5);
        hold on
        
        if ismember('fitTight',P.Results.plot)
            mask = wavetemplatemask(template, Fs, 'noiseThr',0);
            mask(:) = {true(waveLen,1)};
        end
        
        for un = 1:nUnit
            if ismember('fit',P.Results.plot)
                s = s*0;
                s(spkIdx{un}) = 1;
                fh(un) = plot(t,conv(s,template(:,ch,un),'same')/yLim(ch)+yPos(ch),'color',cmap(un,:),'linewidth',1.5);
            else
                spkIdx{un}(spkIdx{un}<waveLen | spkIdx{un}>nDataPoints-waveLen) = [];
                frame = -waveLen/2:waveLen/2-1; %find(mask{ch,un}) - waveLen/2;
                frameMid = frame==0; %round(length(frame)/2);
                if isempty(frame)
                    continue
                end
                txtPos = 1.15*(max(template(:,ch,un))/yLim(ch))+yPos(ch);
                for kk = 1:length(spkIdx{un})
                    fh(un) = plot(t(frame+spkIdx{un}(kk)),template(mask{ch,un},ch,un)/yLim(ch)+yPos(ch),'color',cmap(un,:),'linewidth',2);
                    text(t(frame(frameMid)+spkIdx{un}(kk)),txtPos,num2str(un),'HorizontalAlignment','center','FontSize',10,'color',cmap(un,:))
                end
            end
        end
        
        if ch == nChan
            uNo = num2cell(1:nUnit);
            isPlotted = ~cellfun(@isempty,spkIdx);
            legend(fh(isPlotted),cellfun(@(n) sprintf('unit %i',n),uNo(isPlotted),'uni',false))
        end
    end
%     linkaxes(ax,'x')
    set(gca,'xlim',t([1 end]))
    
    set(gca,'ylim',[-1, yPos(1)+1])
    set(gca,'ytick',flipud(yTickPos))
    set(gca,'yTickLabel',cellfun(@num2str,flipud(num2cell(yTick)),'uni',false))
end

% reconstruction
if ismember('rec',P.Results.plot)
    if holdFig
        figure(figNo(figIdx))
        clf
        figIdx = figIdx+1;
    else
        figure
    end
    hold on
    s = zeros(nDataPoints,1);
    ax = [];
    for ch = 1:nChan
        
        ax = subplot(nChan,1,ch);
        plot(t,double(X(:,ch)),'k','linewidth',1);
        hold on
        
        rec = zeros(nDataPoints,1);
        for un = 1:nUnit
            s = s*0;
            if ~isempty(spkIdx{un})
                s(spkIdx{un}) = 1;
                rec = rec + conv(s,template(:,ch,un),'same');
            end
        end
        plot(t,rec,'c')
    end
    title('original (black) and reconstructed (cyan) signal')
    linkaxes(ax,'x')
    set(gca,'xlim',t([1 end]))
    
%     set(gca,'ylim',[-1, yPos(1)+1])
%     set(gca,'ytick',flipud(yTickPos))
%     set(gca,'yTickLabel',cellfun(@num2str,flipud(num2cell(yTick)),'uni',false))
    
%     for ch = 1:nChan
%         text(gca,1.01*(len+t(end)),yPos(ch),['channel ' num2str(ch)],'fontsize',14)
%     end
end

% residual
if ismember('res',P.Results.plot)
   
   s = zeros(nDataPoints,1);
   res = double(X);
   for un = 1:nUnit
       s = s*0;
       s(spkIdx{un}) = 1;
       for ch = 1:nChan
           res(:,ch) = res(:,ch) - conv(s,template(:,ch,un),'same');
       end
   end
    
   if holdFig
       figure(figNo(figIdx))
       clf
       figIdx = figIdx+1;
   else
       figure
   end
   hold on
   ax = [];
   for ch = 1:nChan
       ax(ch) = subplot(nChan,1,ch);
       plot(t,double(X(:,ch)),'k');
       hold on
       plot(t,res(:,ch),'r');
   end
   title('raw (black) and residual (red)')
   
   linkaxes(ax,'x')
   set(gca,'xlim',t([1 end]))
   
%    set(gca,'ylim',[-1, yPos(1)+1])
%    set(gca,'ytick',flipud(yTickPos))
%    set(gca,'yTickLabel',cellfun(@num2str,flipud(num2cell(yTick)),'uni',false))
   
%    for ch = 1:nChan
%        text(gca,1.01*(len+t(end)),yPos(ch),['channel ' num2str(ch)],'fontsize',14)
%    end
end

% residual energy
if ismember('resen',P.Results.plot)
    
    totResEner = zeros(2,nChan);
    for ch = 1:nChan
        res = X(:,ch);
        for un = 1:nUnit
            s = s*0;
            s(spkIdx{un}) = 1;
            res = res - conv(s,template(:,ch,un),'same');
        end
        totResEner(1,ch) = sum(res.^2);
        totResEner(2,ch) = sum(res(spkRegion(:,ch)>0.5).^2);
    end
    totResEner(1,:) = totResEner(1,:)./sum(X.^2,1);
    for ch = 1:nChan
        totResEner(2,ch) = totResEner(2,ch)./sum(X(spkRegion(:,ch)>0.5,ch).^2);
    end
    totResEner = totResEner';
    
    if ~holdFig
        figure
        bar([resEner{1}, totResEner(:,1)]')
        ylabel('normalized residual energy')
        set(gca,'xtick',1:nUnit+1)
        set(gca,'xticklabel',[cellfun(@(n) sprintf('unit %i',n),num2cell(1:nUnit),'uni',false),{'total'}])
        set(gca,'ylim',[0 max(1,max(max([resEner{1},totResEner(:,1)])))])
        legend(cellfun(@(n) sprintf('chan %i',n),num2cell(1:nChan),'uni',false),'location','northeast')
    end
    
    if holdFig
        figure(figNo(figIdx))
        clf
    else
        figure
    end
    bar([resEner{2}, totResEner(:,2)]')
    ylabel('normalized residual energy (spike regions)')
    set(gca,'xtick',1:nUnit+1)
    set(gca,'xticklabel',[cellfun(@(n) sprintf('unit %i',n),num2cell(1:nUnit),'uni',false),{'total'}])
    set(gca,'ylim',[0 max(1,max(max([resEner{2},totResEner(:,2)])))])
    legend(cellfun(@(n) sprintf('chan %i',n),num2cell(1:nChan),'uni',false),'location','northeast')
end

%% Plot template cross-correlation

% figure
% xcLim = [min(fProj(:)), max(fProj(:))];
% for ii = 1:nUnit
%     for jj = 1:nUnit
%         if ii~=jj
%             subplot(nUnit,nUnit,jj+(ii-1)*nUnit)
%             plot(lags,squeeze(fProj(ii,jj,:)),'k');
%             set(gca,'ylim',xcLim)
%             set(gca,'xlim',lags([1 end]))
%         end
%         if ii == 1 || (ii==2 && jj==1)
%             title(sprintf('f_{%i}',jj))
%         end
%         if mod(jj,nUnit) == 1 || (ii==1 && jj==2)
%             ylabel(sprintf('eta_{%i}',ii))
%         end
%     end
% end

%% Plot raw vs discrimination function (1 plot)

% ax = [];
% figure
% cmap = brewermap(nUnit,'Dark2');
% for ii = 1:nChan
%     ax(ii) = subplot(nChan+1,1,ii);
%     plot(X(:,ii),'k');
% end
% ax(nChan+1) = subplot(nChan+1,1,nChan+1);
% hold on
% for ii = 1:nUnit
%     plot(d{ii},'color',cmap(ii,:))
% end
% yline(thr,'color',0.8*[1 1 1]);
% linkaxes(ax,'x')

%% Plot raw vs discrimination function (multi-plot)

% cmap = brewermap(nUnit,'Dark2');
% for ii = 1:nChan
%     ax = [];
%     figure
%     ax(1) = subplot(211);
%     plot(X(:,ii),'k');
%     title(sprintf('channel %i',ii))
%     ax(2) = subplot(212);
%     hold on
%     for jj = 1:nUnit
%         plot(d{jj},'color',cmap(jj,:))
%     end
%     yline(thr,'color',0.8*[1 1 1]);
%     linkaxes(ax,'x')
% end
        