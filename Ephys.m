%% Electrophysiologyiological Data Class
%   For storing, viewing, and processing electrophysiological data
%
%   CONSTRUCTOR
%       E = Ephys(Fs,data,alignIndex,chanNo);
%
%   PROPERTIES
%       Fs <scalar> - data sample rate in Hz (default: 1000)
%       data <numeric array> - timeseries data of dimensionality <= 3
%       alignIdx <integer> - alignment index (default: 1)
%       chanNo <
%
%   METHODS
%   ~type 'Ephys.man' for full user-manual
classdef Ephys
    properties (Access = public)
        Fs = 1e3            % sample frequency (Hz)
        data                % timeseries data (samples x channels x trials)
        alignIndex = 1      % aligment index (scalar integer)
        chanNo              % channel number (1D integer array)
    end
    properties (SetAccess = private)
        time                % time vector
        nSamples            % number of data points
        nChannels           % number of channels
        nTrials             % number of trials
        trialNo             % trial number
    end
    methods (Static)
        function obj = init_from_nsx(filePath)
            % initializes object from NSx file
            nsx = openNSx(filePath);
            obj = Ephys;
            obj.Fs = nsx.MetaTags.SamplingFreq;
            if iscell(nsx.Data)
                nsx.Data = cat(2,nsx.Data{:});
            end
            obj.data = nsx.Data(1:nsx.MetaTags.ChannelCount-1,:)';
        end
    end
    methods
        % -----------------------------------------------------------------
        % CONSTRUCTOR
        % -----------------------------------------------------------------
        function obj = Ephys(Fs, data, alignIndex, chanNo)
            if nargin > 0 && isnumeric(Fs) && isscalar(Fs)
                obj.Fs = Fs;
            end
            if nargin > 1 && isnumeric(data) && ndims(data)<=3
                    obj.data = data;
            end
            if nargin > 2 && isnumeric(alignIndex) && isscalar(alignIndex) && alignIndex>0
                obj.alignIndex = alignIndex;
            end
            if nargin > 3 && isnumeric(chanNo)
                obj.chanNo = chanNo;
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
        function obj = set.data(obj, data)
            if isnumeric(data) && ndims(data)<=3
                obj.data = data;
            else
                error('Invalid data assignment')
            end
            obj = obj.update_time;
            obj = obj.update_trial;
        end
        function obj = set.alignIndex(obj, alignIndex)
            if isscalar(alignIndex) && isnumeric(alignIndex)
                obj.alignIndex = alignIndex;
            else
                error('Invalid alignment index assignment')
            end
            obj = obj.update_time;
        end
        function obj = set.chanNo(obj, chanNo)
            if isnumeric(chanNo)
                obj.chanNo = chanNo;
            else
                error('Invalid channel label assignment')
            end
        end
        % -----------------------------------------------------------------
        % GET METHODS
        % -----------------------------------------------------------------
        function nPts = get.nSamples(obj)
            nPts = size(obj.data,1);
        end
        function nChan = get.nChannels(obj)
            nChan = size(obj.data,2);
        end
        function nTr = get.nTrials(obj)
            nTr = size(obj.data,3);
        end
        function chNo = get.chanNo(obj)
            if isempty(obj.chanNo)
                chNo = 1:obj.nChannels;
            else
                chNo = obj.chanNo;
            end
        end
        % -----------------------------------------------------------------
        % SPECIAL
        % -----------------------------------------------------------------
        function disp(obj)
            txt = '%i-channel ephys data sampled at %.3f kHz for %.3f %s\n';
            durSec = obj.nSamples/obj.Fs;
            if durSec < 60
                dur = durSec;
                durUnit = 'sec';
                
            elseif durSec > 60 && durSec < 3600
                dur = durSec/60;
                durUnit = 'min';
                
            else
                dur = durSec/3600;
                durUnit = 'hrs';
            end
           fprintf(txt,obj.nChannels,obj.Fs/1e3,dur,durUnit);
        end
        function Eph = plus(E1,E2)
            assert(E1.Fs == E2.Fs)
            assert(size(E1,2) == size(E2,2))
            Eph = Ephys(E1.Fs);
            Eph.data = cat(1, E1.data, E2.data);
            Eph.chanNo = [E1.chanNo, E2.chanNo];
        end
        % -----------------------------------------------------------------
        % DATA PARSING
        % -----------------------------------------------------------------
        function [obj,tIdx] = range(obj,tLim)
            % truncates data between time limits
            assert(tLim(1)>=0 && tLim(2)<=obj.time(end), 'Time range out of bounds')
            tIdx = (obj.time>=tLim(1) & obj.time<=tLim(2));
            obj.data = obj.data(tIdx,:,:);
            obj.alignIndex = 1-find(tIdx,1);
        end
        function obj = chan(obj,idx)
            % pulls channel indices
            assert(all(idx<=obj.nChannels),'Channel out of bounds')
            obj.chanNo = obj.chanNo(idx);
            obj.data = obj.data(:,idx,:);
        end
        function obj = trial(obj,num)
            % selects specified trials from data
            if isscalar(num)
                assert(num>=1 && num<=obj.nTrials)
            elseif ischar(num) && strcmp(num,'rand')
                num = randi(obj.nTrials);
            else
                error('Provide a scalar trial number or specify as ''rand''.')
            end
            obj.data = obj.data(:,:,num);
            obj.trialNo = num;
        end
        function [obj,v] = cat_trials(obj,num)
            % concatenates trials
            if nargin == 1
                num = 1:obj.nTrials;
            end
            v = cellfun(@(x) obj.data(:,:,x), num2cell(num(:)),'uni',false);
            v = cell2mat(v);
            obj.data = v;
            obj.trialNo = num;
        end
        % -----------------------------------------------------------------
        % PROCESSING
        % -----------------------------------------------------------------
        function obj = rect(obj)
            % rectifies data
            for ii = 1:length(obj)
                obj(ii).data = abs(obj(ii).data);
            end
        end
        function obj = filt(obj,type,arg1,arg2)
            % filter data
            % .filt(type,n,w) if type is 'bandpass','low', or 'high' applies
            %   a corresponding nth-order butterworth filter with cutoffs w [Hz]
            % .filt('diff',n) applies an nth order differentiating filter
            % .filt('gau',sd) smooths with a Gaussian filter with standard
            %   deviation sd [seconds]
            % .filt('box',wid) smooths with a boxcar filter of width wid [seconds]
            for ii = 1:length(obj)
                switch type
                    case {'bandpass','low','high'}
                        n = arg1;
                        w = arg2;
                        
                        [b,a] = butter(n, w/(obj(ii).Fs/2), type);
                        dataClass = class(obj(ii).data);
                        
                        xF = mat2cell(obj(ii).data, obj(ii).nSamples,...
                            ones(1,obj(ii).nChannels,1), ones(1,1,obj(ii).nTrials));
                        xF = cellfun(@(y) filtfilt(b,a,double(y)), xF, 'uni', false);
                        xF = cell2mat(xF);
                        xF = eval(sprintf('%s(xF);',dataClass));
                        obj(ii).data  = xF;
                        
                    case 'diff'
                        n = arg1;
                        obj(ii).data = diff(obj(ii).data,n,1);
                        
                    case 'gau'
                        sd = arg1;
                        obj(ii).data = smooth1D(obj(ii).data,obj(ii).Fs,'gau','sd',sd,'dim',1);
                        
                    case 'box'
                        wid = arg1;
                        obj(ii).data = smooth1D(obj(ii).data,obj(ii).Fs,'box','widBox',wid,'dim',1);
                        
                    otherwise
                        error('Unrecognized filter type')
                end
            end
        end
        function [obj,sigma] = normalize(obj)
            % normalize data
            % [obj,sigma] = obj.norm normalizes data by the estimated
            %   standard deviation of the noise, sigma
            sigma = cell(length(obj),1);
            for ii = 1:length(obj)
                sigma{ii} = median(abs(double(obj(ii).data)),1)/0.6745;
                obj(ii).data = double(obj(ii).data)./sigma{ii};
            end
            if length(obj) == 1
                sigma = sigma{1};
            end
        end
        function [Cinv, Cn] = noise_cov(obj,t0,waveDur)
            xNoise = EMG.range(t0+[0,1-1/obj.Fs]).data;
            Cn = noisecov(obj.data,obj.Fs,waveDur,'xNoise',xNoise);
            Cinv = Cn^(-1);
        end
        function [trMean, trSte, trStd] = trial_avg(obj)
            % trial average
            trMean = mean(double(obj.data),3);
            trStd = std(double(obj.data),[],3);
            trSte = trStd/sqrt(obj.nTrials-1);
            
            % convert to native data class
            dataClass = class(obj.data);
            if ~isa(obj.data,'double')
                eval(sprintf('trMean = %s(trMean);',dataClass))
                eval(sprintf('trSte = %s(trSte);',dataClass))
            end
        end
        function [res,vEst] = residual(obj,Neu)
            % compute residual
            assert(isa(Neu,'Neuron'))
            vEst = zeros(obj.nSamples, obj.nChannels);
            idxLim = [1 length(obj.time)];
            for ch = 1:obj.nChannels
                for un = 1:Neu.nUnits
                    si = find(Neu.spikes(:,un)); % + obj.alignIndex;
                    si = si(si>idxLim(1) & si<idxLim(2));
                    x = zeros(obj.nSamples,1);
                    x(si) = 1;
                    vEst(:,ch) = vEst(:,ch) + conv(x,Neu.waveform(:,ch,un),'same');
                end
            end
            eval(sprintf('yEst = %s(yEst);',class(obj.data)))
            res = obj.data - vEst;
        end
        % -----------------------------------------------------------------
        % WRAPPER FUNCTIONS
        % -----------------------------------------------------------------
        function Spk = detect_spikes(obj,varargin)
            % detect spike events and return a Spike object
            isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;
            P = inputParser;
            addOptional(P, 'thresh', 4, @(x) isscalarnum(x,-Inf,Inf))
            addParameter(P, 'sym', true, @islogical)
            addParameter(P, 'sigma', [], @(x) isempty(x) || isnumeric(x))
            addParameter(P, 'minWid', 1e-3, @(x) isscalarnum(x,0,1))
            addParameter(P, 'filt', true, @islogical)
            addParameter(P, 'normalized', false, @islogical)
            parse(P,varargin{:})
            params = [P.Parameters(:),struct2cell(P.Results)];
            params = reshape(params',1,numel(params));
            Spk = repmat(Spike(obj.Fs),obj.nChannels,1);
            idx = findpulses(double(obj.data),obj.Fs,params{:});
            if ~iscell(idx)
                idx = {idx};
            end
            for ii = 1:obj.nChannels
                Spk(ii).index = idx{ii};
            end
        end
        % -----------------------------------------------------------------
        % PLOTTING
        % -----------------------------------------------------------------
        function fh = set_fig(~,code)
            % set figure
            if isempty(code)
                return
            end
            switch code
                case 'clf'
                    clf
                    fh = gcf;
                case 'new'
                    fh = figure;
            end
        end
        function plot(obj,varargin)
            % plot
            P = inputParser;
            addParameter(P,'fig',[],@(x) ischar(x) && ismember(x,{'clf','new'}))
            addParameter(P,'Spk',[],@(x) isa(x,'Spike'))
            addParameter(P,'Neu',[],@(x) isa(x,'Neuron'))
            addParameter(P,'sigma',[],@isnumeric)
            addParameter(P,'txt',[],@ischar)
            parse(P,varargin{:})
            
            obj.set_fig(P.Results.fig);
            
            ax = [];
            for ii = 1:obj.nChannels
                if obj.nChannels > 1
                    ax(ii) = subplot(obj.nChannels,1,ii);
                    if ~isempty(P.Results.sigma)
                        hold on
                        yline(P.Results.sigma,'r');
                        yline(-P.Results.sigma,'r');
                    end
                end
                plot(obj.time,obj.data(:,ii),'k')
                title(sprintf('channel %i',obj.chanNo(ii)))
            end
            xlabel('time (s)')
            if obj.nChannels > 1
                linkaxes(ax,'x')
            end
            set(gca,'xlim',obj.time([1 end]))
            if ~isempty(P.Results.txt)
                th = text(0,0,P.Results.txt);
                th.Units = 'normalized';
                th.Position = [-0.05,-0.1];
                th.FontSize = 15;
                th.BackgroundColor = 'w';
            end                
            
            % overlay spike times
            Spk = P.Results.Spk;
            if ~isempty(Spk)
                spkLim = round(obj.Fs * obj.time([1 end]));
                for ii = 1:obj.nChannels
                    if obj.nChannels > 1
                        subplot(obj.nChannels,1,ii)
                    end
                    hold on
                    spkIdx = Spk(ii).index(Spk(ii).index>=spkLim(1) & Spk(ii).index<=spkLim(2));
                    spkIdx = 1+spkIdx-spkLim(1);
                    spkAmp = obj.data(spkIdx,ii);
                    plot(obj.time(spkIdx),1.25*spkAmp,'r*')
                end
            end
            
            % overlay neuronal waveforms
            Neu = P.Results.Neu;
            if ~isempty(Neu)
                M = plotmarkers();
                win = (0:size(Neu.waveform,1)-1) - size(Neu.waveform,1)/2;
                idxLim = [1+length(win), length(obj.time)-length(win)];
                wEner = squeeze(sum(Neu.waveform.^2,1));
                
                for ii = 1:obj.nChannels
                    if obj.nChannels > 1
                        subplot(obj.nChannels,1,ii)
                    end
                    hold on
                    fh = zeros(1,Neu.nUnits);
                    for un = 1:Neu.nUnits
                        if wEner(ii,un) == 0
                            continue
                        end
                        si = find(Neu.spikes(:,un)); % + obj.alignIndex;
                        si = si(si>idxLim(1) & si<idxLim(2));
                        for kk = 1:length(si)
                            fh(un) = plot(obj.time(win+si(kk)),Neu.waveform(:,ii,un),M.line{un+1});
                        end
                    end
                    if ii == 1
                        unitNames = cellfun(@(n) sprintf('unit %i',n), num2cell(1:Neu.nUnits), 'uni', false);
                        legend(fh(fh>0), unitNames(fh>0))
                    end
                end
            end          
        end
    end
    methods (Access = private)
        function obj = update_time(obj)
            obj.time = ((1:obj.nSamples)'-obj.alignIndex)/obj.Fs;
        end
        function obj = update_trial(obj)
            if isempty(obj.trialNo)
                obj.trialNo = 1:size(obj.data,3);
            end
        end
    end            
end