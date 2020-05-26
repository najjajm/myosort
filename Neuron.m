%% Neuronal Data Class
classdef Neuron
    properties
        Fs = 1e3
        spikes
        alignIndex = 1
        waveform
    end
    properties (SetAccess = private)
        rate
        psth
        nDataPoints
        nUnits
        nTrials
    end
    properties (Constant, Hidden)
        filtSamplingFrequency = 1e3
        gaussSD = 25e-3
    end
    methods (Static)
        % -----------------------------------------------------------------
        % INITIALIZE FROM SPIKE OBJECT
        % -----------------------------------------------------------------
        function obj = init_from_spike(Spk)
            Fs = Spk(1).Fs;
            waveLen = Spk(1).waveLength;
            obj = Neuron(Fs);
            % read spike object
            nChan = length(Spk);
            unitNo = arrayfun(@(x) unique(x.label(x.label>0)), Spk, 'uni', false);
            unqUnitNo = unique(cell2mat(unitNo));
            nTotUnits = length(unqUnitNo);
            obj.waveform = zeros(waveLen,nChan,nTotUnits);
            % copy waveforms
            for un = 1:nTotUnits
                for ch = 1:nChan
                    w = Spk(ch).waveform(Spk(ch).label==unqUnitNo(un),:);
                    if ~isempty(w)
                       obj.waveform(:,ch,un) = mean(double(w))'; 
                    end
                end
            end
            % copy spike data
            s = cell(1,nTotUnits);
            for un = 1:nTotUnits
               ch = find(cellfun(@(x) ismember(unqUnitNo(un),x), unitNo));
               if length(ch) == 1
                   s{un} = sparse(Spk(ch).index(Spk(ch).label==unqUnitNo(un)),1,true,Spk(1).nDataPoints,1);
               else
                  sTmp = zeros(Spk(1).nDataPoints,length(ch));
                  for ii = 1:length(ch)
                      sTmp(Spk(ch(ii)).index(Spk(ch(ii)).label==unqUnitNo(un)),ii) = 1;
                  end
                  sF = smooth1D(sum(sTmp,2), Fs, 'gau', true, 'sd', waveLen/(2*Fs));
                  [~,pkLoc] = findpeaks(sF, 'MinPeakHeight', 0.9);
                  s{un} = sparse(pkLoc,1,true,Spk(1).nDataPoints,1);
               end
            end
            obj.spikes = cell2mat(s);
        end
    end
    methods
        % -----------------------------------------------------------------
        % CONSTRUCTOR
        % -----------------------------------------------------------------
        function obj = Neuron(Fs, spikes, alignIndex)
            
            if nargin > 0 && isnumeric(Fs) && isscalar(Fs)
                obj.Fs = Fs;
            end
            if nargin > 1 && (isa(spikes,'ndSparse')||issparse(spikes)) && ndims(spikes)<=3
                obj.spikes = spikes;
            end
            if nargin > 2 && isnumeric(alignIndex) && isscalar(alignIndex) && alignIndex>0
                obj.alignIndex = alignIndex;
            end
        end
        
        % -----------------------------------------------------------------
        % SET METHODS
        % -----------------------------------------------------------------
        function obj = set.Fs(obj, Fs)
            if isnumeric(Fs) && isscalar(Fs)
                obj.Fs = Fs;
            else
                error('Invalid sampling frequency assignment')
            end
        end
        function obj = set.spikes(obj, spikes)
            if (isa(spikes,'ndSparse')||issparse(spikes)) && ndims(spikes)<=3
                obj.spikes = spikes;
                obj = obj.update_rate;
                obj = obj.update_psth;
            else
                if ~isempty(spikes)
                    error('Invalid spike array assignment')
                else
                    obj.spikes = sparse([]);
                end
            end
        end
        function obj = set.alignIndex(obj, alignIndex)
            if isnumeric(alignIndex) && isscalar(alignIndex) && alignIndex>0
                obj.alignIndex = alignIndex;
            else
                error('Invalid zero index assignment')
            end
        end
        
        % -----------------------------------------------------------------
        % GET METHODS
        % -----------------------------------------------------------------
        function nPts = get.nDataPoints(obj)
            nPts = size(obj.spikes,1);
        end
        function nUnit = get.nUnits(obj)
            nUnit = size(obj.spikes,2);
        end
        function nTr = get.nTrials(obj)
            nTr = size(obj.spikes,3);
        end
        
        % -----------------------------------------------------------------
        % SELECT TRIAL
        % -----------------------------------------------------------------
        function obj = trial(obj,trNo)
            assert(trNo>=1 && trNo<=obj.nTrials)
            obj.spikes = obj.spikes(:,:,trNo);
        end
        
        % -----------------------------------------------------------------
        % DISPLAY
        % -----------------------------------------------------------------
        
        % -----------------------------------------------------------------
        % DISPLAY
        % -----------------------------------------------------------------
        
        % -----------------------------------------------------------------
        % PROCESSING
        % -----------------------------------------------------------------
        
        % trial time vector
        function t = time(obj,len,Fs)
            if nargin < 2
                len = obj.nDataPoints;
            end
            if nargin < 3
                Fs = obj.Fs;
            end
            alignIdx = round(obj.alignIndex*(Fs/obj.Fs));
            t = ((1:len)-alignIdx)/Fs;
        end
        
        % normalize
        function obj = normalize(obj,varargin)
            P = inputParser;
            addParameter(P,'soft',false,@islogical)
            addParameter(P,'normFact',1.1,@(x) isnumeric(x) && x>=1)
            parse(P,varargin{:})
            
            rMax = max(obj.psth,[],1);
            if P.Results.soft
                rMax = rMax * P.Results.normFact;
            end
            obj.psth = obj.psth./rMax;           
        end
        
        % -----------------------------------------------------------------
        % PLOTTING
        % -----------------------------------------------------------------
        
        % spike rasters
        function plot_rast(obj,varargin)
            P = inputParser;
            addParameter(P,'unit',1:obj.nUnits,@(x) isnumeric(x) && min(x)>0 && max(x)<=obj.nUnits)
            addParameter(P,'buffer',0.5,@(x) isnumeric(x) && x>0)
            addParameter(P,'legend',true,@islogical)
            addParameter(P,'tickDur',[],@(x) isempty(x) || isnumeric(x))
            addParameter(P,'lineWidth',1,@isnumeric)
            parse(P,varargin{:})
            
            s = obj.spikes;
            t = obj.time;
            
            cmap = obj.unit_colors(-0.2);
            hold on
            
            unitNo = P.Results.unit;
            nUnit = length(unitNo);
            unitNo = fliplr(unitNo);
            
            [trialNo,lineNo] = deal(repmat((1:obj.nTrials)',1,nUnit)); 
            
            ySpace = ceil(P.Results.buffer*obj.nTrials);
            lineBuff = (obj.nTrials+ySpace)*(0:nUnit-1);
            lineBuff(2:end) = lineBuff(2:end)-1;
            lineNo = lineNo + lineBuff;
            
            trialNo = trialNo(:);
            lineNo = fliplr(lineNo);
            lineNo = lineNo(:);
            
            ySampIdx = repmat(round(linspace(1,obj.nTrials,min(obj.nTrials,3)))',1,nUnit);
            ySampIdx = ySampIdx + (obj.nTrials)*(0:nUnit-1);
            ySampIdx = fliplr(ySampIdx);
            ySampIdx = ySampIdx(:);
            
            if ~isempty(P.Results.tickDur)
               tickLen = round(P.Results.tickDur*obj.Fs); 
            end
            
            fh = zeros(1,nUnit);
            for ii = 1:nUnit
                for jj = 1:obj.nTrials
                    sIdx = find(s(:,unitNo(ii),jj));
                    if isempty(sIdx)
                        continue
                    end
                    lNo = lineNo(jj+obj.nTrials*(ii-1));
                    if isempty(P.Results.tickDur)
                        fh(ii) = plot(t(sIdx),lNo*ones(1,length(sIdx)),'.','color',cmap(unitNo(ii),:));
                    else
                        sIdx(sIdx+(tickLen-1) > length(t)) = [];
                        for kk = 1:length(sIdx)
                            fh(ii) = plot(t(sIdx(kk)+[0 tickLen-1]),lNo*[1 1],'color',cmap(unitNo(ii),:),'linewidth',P.Results.lineWidth);
                        end
                    end
                end
            end
            set(gca,'xlim',t([1 end]))
            set(gca,'ylim',[0 1+max(lineNo)])
            
            set(gca,'ytick',lineNo(ySampIdx))
            set(gca,'yticklabel',cellfun(@num2str,num2cell(trialNo(ySampIdx)),'uni',false))
            ylabel('trial')
            
            if P.Results.legend
                legTxt = cellfun(@(n) sprintf('MU %i',n), num2cell(unitNo), 'uni',false);
                legend(fh(fh~=0),legTxt(fh~=0),'location','best')
            end
        end
        
        % firing rates
        function plot_rates(obj,varargin)
            P = inputParser;
            addParameter(P,'unit',1:obj.nUnits,@(x) isnumeric(x) && min(x)>0 && max(x)<=obj.nUnits)
            addParameter(P,'trial',1:obj.nTrials,@(x) isnumeric(x) && min(x)>0 && max(x)<=obj.nTrials)
            parse(P,varargin{:})
            
            fr = obj.rate;
            t = obj.time(size(fr,1),obj.filtSamplingFrequency);
            
            nUnit = length(P.Results.unit);
            nCol = ceil(sqrt(nUnit));
            nRow = ceil(nUnit/nCol);
            nTrial = length(P.Results.trial);
            shade = linspace(0.25,0.75,nTrial);
            for ii = 1:nUnit
                if nUnit > 1
                    subplot(nRow,nCol,ii)
                end
                hold on
                for jj = 1:nTrial
                    plot(t,fr(:,P.Results.unit(ii),P.Results.trial(jj)),'color',shade(jj)*ones(1,3))
                end
            end
            set(gca,'xlim',t([1 end]))
            xlabel('time (s)')
            ylabel('firing rate (Hz)')
        end
        
        % peri-stimulus time histograms
        function plot_psth(obj,varargin)
            P = inputParser;
            addParameter(P,'unit',1:obj(1).nUnits,@(x) isnumeric(x) && min(x)>0 && max(x)<=obj(1).nUnits)
            addParameter(P,'cond',1:length(obj),@isnumeric)
            addParameter(P,'fig',[],@(x) ischar(x) && ismember(x,{'clf','new'}))
            addParameter(P,'label',[],@ischar)
            addParameter(P,'tLim',[-Inf,Inf],@(x) isnumeric(x) && length(x)==2)
            addParameter(P,'legend',true,@islogical)
            parse(P,varargin{:})
            
            if ~isempty(P.Results.fig)
                if strcmp(P.Results.fig,'clf')
                    clf
                elseif strcmp(P.Results.fig,'new')
                    figure
                end
            end
            
            nUnit = size(obj(1).psth,2);   
            cmap = obj.unit_colors;
            nCond = length(P.Results.cond);
            nCol = ceil(sqrt(nCond));
            nRow = ceil(nCond/nCol);
            
            absRMax = max(arrayfun(@(x) max(max(x.psth(:,:,1))),obj));
            for ii = 1:nCond
                if nCond > 1
                    subplot(nRow,nCol,ii)
                end
                cNo = P.Results.cond(ii);
                
                y = obj(cNo).psth;
                t = obj(cNo).time(size(y,1),obj(cNo).filtSamplingFrequency);
                
                tIdx = t>=P.Results.tLim(1) & t<=P.Results.tLim(2);
                t = t(tIdx);
                y = y(tIdx,:,:);
                
                fh = zeros(1,nUnit);
                hold on
                for jj = 1:length(P.Results.unit)
                    unitNo = P.Results.unit(jj);
                    fh(jj) = plot(t,y(:,unitNo,1),'color',cmap(unitNo,:),'linewidth',3);
                    patch([t,fliplr(t)], [permute(sum(y(:,unitNo,:),3),[2 1 3]),...
                        fliplr(permute(-diff(y(:,unitNo,:),[],3),[2 1 3]))], cmap(unitNo,:),...
                        'EdgeAlpha',0, 'FaceAlpha',0.125)
                end
                set(gca,'xlim',t([1 end]))
                set(gca,'ylim',[0 max(1,round(1.1*absRMax))]) % max(1,1.1*max(max(sum(y,3))))])
                xlabel('time (s)')
                ylabel('firing rate (Hz)')
                title(sprintf('condition %i',cNo))
                
                if ii == 1 && P.Results.legend
                    legText = cellfun(@(n) sprintf('unit %i',n),num2cell(P.Results.unit),'uni',false);
                    legend(fh(fh>0),legText(fh>0),'location','NorthWest')
                end
            end
            if ~isempty(P.Results.label)
                subplot(nRow,nCol,1+(nRow-1)*nCol)
                th = text(0,0,P.Results.label);
                th.Units = 'normalized';
                th.Position = [-.1 -.2];
                th.FontSize = 12;
            end
        end
    end
    methods (Access = private)
        % -----------------------------------------------------------------
        % MOTOR UNIT COLOR MAP
        % -----------------------------------------------------------------
        function cmap = unit_colors(obj,brightness)   
            if nargin == 1
                brightness = -0.1;
            end
            cmap = brewermap(obj(1).nUnits,'Spectral');
            cmap = max([0 0 0],cmap + brightness); 
        end
        % -----------------------------------------------------------------
        % UPDATE FIRING RATES
        % -----------------------------------------------------------------
        function obj = update_rate(obj)
            Fs = obj.filtSamplingFrequency;
            spk = mat2cell(obj.spikes,obj.nDataPoints,ones(1,obj.nUnits,1),ones(1,1,obj.nTrials));
            % downsample spikes
            if obj.Fs ~= Fs
                t = obj.time();
                tNew = t(1):1/Fs:t(end);
                tSpk = cellfun(@(s) t(s), spk, 'uni', false);
                spk = cellfun(@(tS) hist(tS,tNew), tSpk, 'uni', false);
                spk = cellfun(@(x) x(:), spk, 'uni', false);
            end
            % filter
            fr = cell2mat(cellfun(@(s) smooth1D(double(full(s)),Fs,...
                'gau','sd',obj.gaussSD), spk, 'uni', false));
            obj.rate = Fs * fr;
        end
        % -----------------------------------------------------------------
        % UPDATE PSTHs
        % -----------------------------------------------------------------
        function obj = update_psth(obj)
            obj.psth = cat(3,mean(obj.rate,3),std(obj.rate,[],3)/sqrt(obj.nTrials));
        end
    end
end