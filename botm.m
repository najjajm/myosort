%% BOTM Bayes optimal template matching
% Function details
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
addParameter(P, 'sic', true, @islogical)
addParameter(P, 'plot', {''}, @(x) iscell(x) && all(ismember(x,{'','fit','res','resen','rec','dis'})))
addParameter(P, 'cmap', 'Spectral', @ischar)
addParameter(P, 'verbose', false, @islogical)

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
lags = (1-waveLen):(waveLen-1);

% whitened filter energy
filtEner = zeros(nChan,nUnit);
for ii = 1:nChan
    filtEner(ii,:) = cellfun(@(a,b) -1/2*a(:,ii)'*b(:,ii), w, u);
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
    XFT = fft(X(:,ch),nPad,1);
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

if P.Results.verbose
    fprintf('Runtime: %.2f min\n',toc(t0)/60)
end

%% Subtractive interference cancellation

if P.Results.verbose
    fprintf('Running SIC. ')
    t0 = tic;
end

if P.Results.sic && nnz(~cellfun(@isempty,spkIdx)) > 1
    
    % extract the discriminant at each spike
    dSpk = (thr-eps) * ones(nDataPoints,nUnit);
    for un = 1:nUnit
        dSpk(spkIdx{un},un) = d(spkIdx{un},un);
    end
    
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
                if dSpk(t1,u1) < thr
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
                    dSpk(t2,u2) = dSpk(t2,u2) - fProj(u2,u1,lags==tau);
                end
            end
        end
        
        % update spike estimate
        spkIdx = cell(1,nUnit);
        for un = 1:nUnit
            [~,spkIdx{un}] = findpeaks(dSpk(:,un),'MinPeakHeight',thr);
            spkIdx{un} = spkIdx{un}(spkIdx{un}>waveLen & spkIdx{un}<nDataPoints-waveLen);
        end
    end
end

if P.Results.verbose
    fprintf('Runtime: %.2f min\n',toc(t0)/60)
end

%% Residual energy

if P.Results.verbose
    fprintf('Computing residual energy. ')
    t0 = tic;
end

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

%% Plotting

cmap = brewermap(nUnit,P.Results.cmap);
% cmap = flipud(jet);
% cmap = cmap(round(linspace(1,size(cmap,1),nUnit)),:);
cmap = max([0 0 0],cmap-.2);

t = (1:nDataPoints)/Fs;

% discriminant
if ismember('dis',P.Results.plot)
    figure
    ax = [0 0];
    ax(1) = subplot(211);
    hold on
    for ii = 1:nUnit
        plot(d(:,ii),'color',cmap(ii,:));
    end
    set(gca,'xlim',[1 size(d,1)])
    plot(get(gca,'xlim'),thr*[1 1],'color',0.8*[1 1 1])
    title(sprintf('Fisher''s Discriminant\nraw'))
    legend(cellfun(@(n) sprintf('unit %i',n),num2cell(1:nUnit),'uni',false))
    ax(2) = subplot(212);
    hold on
    for ii = 1:nUnit
        plot(dSpk(:,ii),'color',cmap(ii,:),'linewidth',1.5)
    end
    title('post-SIC')
    set(gca,'xlim',[1 size(d,1)])
    plot(get(gca,'xlim'),thr*[1 1],'color',0.8*[1 1 1])
    linkaxes(ax,'x')
end

% single unit fits
if ismember('fit',P.Results.plot)    
    figure
    ax = zeros(1,nChan);
    s = zeros(nDataPoints,1);
    for ch = 1:nChan
        
        ax(ch) = subplot(nChan,1,ch);
        fh = plot(t,X(:,ch),'k');
        fh.Color(4) = 0.25;
        hold on
        
        fh = [];
        len = waveLen/2;
        frame = -len/2:len/2-1;
        for un = 1:nUnit
%             for kk = 1:length(spkIdx{un})
%                 fh(un) = plot(t(frame+spkIdx{un}(kk)),template(frame+waveLen/2,ch,un),'color',cmap(un,:),'linewidth',1.5);
%             end
            s = s*0;
            s(spkIdx{un}) = 1;
            fh(un) = plot(t,conv(s,template(:,ch,un),'same'),'color',cmap(un,:));
        end
        
        ylabel(sprintf('channel %i',ch))
        if ch == 1
            title('single unit fits')
        end
        if ch == nChan
            uNo = num2cell(1:nUnit);
            isPlotted = ~cellfun(@isempty,spkIdx);
            legend(fh(isPlotted),cellfun(@(n) sprintf('unit %i',n),uNo(isPlotted),'uni',false))
        end
    end
    linkaxes(ax,'x')
    set(gca,'xlim',t([1 end]))
end

% reconstruction
if ismember('rec',P.Results.plot)
    figure
    ax = zeros(1,nChan);
    s = zeros(nDataPoints,1);
    for ch = 1:nChan
        
        ax(ch) = subplot(nChan,1,ch);
        plot(X(:,ch),'k','linewidth',1);
        hold on
        
        rec = zeros(nDataPoints,1);
        for un = 1:nUnit
            s = s*0;
            if ~isempty(spkIdx{un})
                s(spkIdx{un}) = 1;
                rec = rec + conv(s,template(:,ch,un),'same');
            end
        end
        plot(rec,'c')
        
        if ch == 1
            title('original (black) and reconstructed (cyan) signal')
        end
        ylabel(sprintf('channel %i',ch))
    end
    linkaxes(ax,'x')
%     sgtitle('original (black) and reconstructed (cyan) signal')
    set(gca,'xlim',[1 nDataPoints])
end

% residual
if ismember('res',P.Results.plot)
   
   s = zeros(nDataPoints,1);
   res = X;
   for un = 1:nUnit
       s = s*0;
       s(spkIdx{un}) = 1;
       for ch = 1:nChan
           res(:,ch) = res(:,ch) - conv(s,template(:,ch,un),'same');
       end
   end
    
   figure
   ax = zeros(1,nChan);
   for ch = 1:nChan
       ax(ch) = subplot(nChan,1,ch);
       plot(t,X(:,ch),'k');
       hold on
       plot(t,res(:,ch),'r');
       
       if ch == 1
           title('raw (black) and residual (red)')
       end
       ylabel(sprintf('channel %i',ch))
   end
   linkaxes(ax,'x')
   set(gca,'xlim',t([1 end]))
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
    
    figure
    bar([resEner{1}, totResEner(:,1)]')
    ylabel('normalized residual energy')
    set(gca,'xtick',1:nUnit+1)
    set(gca,'xticklabel',[cellfun(@(n) sprintf('unit %i',n),num2cell(1:nUnit),'uni',false),{'total'}])
    set(gca,'ylim',[0 max(1,max(max([resEner{1},totResEner(:,1)])))])
    legend(cellfun(@(n) sprintf('chan %i',n),num2cell(1:nChan),'uni',false),'location','northeast')
    
    figure
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
%         plot(dSpk{jj},'color',cmap(jj,:))
%     end
%     yline(thr,'color',0.8*[1 1 1]);
%     linkaxes(ax,'x')
% end
        