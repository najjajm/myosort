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
                error('Invalid spike array assignment')
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
        
        % -----------------------------------------------------------------
        % PLOTTING
        % -----------------------------------------------------------------
        
        % action potential waveform
        function plot_waveform(obj)
            [waveLen,nChan,nUnit] = size(obj.waveform);
            tW = (1:waveLen)/obj.Fs;
            waveMin = permute(min(obj.waveform,[],1), [2 3 1]);
            waveMax = permute(max(obj.waveform,[],1), [2 3 1]);
            waveRange = diff(cat(3,waveMin,waveMax),[],3)';
            yBound = max(waveRange,[],2);
            yOffset = cumsum(yBound);
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
            
            nUnit = size(obj(1).waveform,3);
%             cmap = flipud(jet);
%             cmap = cmap(round(linspace(1,size(jet,1),nUnit)),:);
%             
%             % temporary (darker yellow)
%             cmap(6,:) = [252 238 33]/255;
            
            cmap = brewermap(nUnit,'Spectral');
            cmap = max([0 0 0],cmap-.1);
            
            nCond = length(P.Results.cond);
            nCol = ceil(sqrt(nCond));
            nRow = ceil(nCond/nCol);
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
                set(gca,'ylim',[0 max(1,1.1*max(max(sum(y,3))))])
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