%% Electrophysiologyiological Data Class
%
% CONSTRUCTOR
%   Eph = Ephys(Fs, data, alignIndex, chanNo)
%
% PROPERTIES (public access)
%   Fs <scalar>: sampling frequency in Hz
%   data <numeric>: array of dimensions (samples x channels x trials)
%   alignIndex <scalar>: alignment point (i.e. sample corresponding to t=0)
%       (default: 1)
%   chanNo <numeric>: vector of channel numbers (default: 1:size(data,2))
%
% PROPERTIES (private set access)
%   time: vector in seconds
%   nSamples: number of data points
%   nChannels: channel count
%   nTrials: trial count
%   trialNo: trial number
%
% METHODS (static)
%   Eph = Ephys.init_from_nsx(filePath) initalizes an Ephys object by 
%       reading an NSx file from specified path
%
% METHODS (public access)
%   plus(E1,E2) overloads the normal "plus" operation to define adding 
%       two Ephys objects (i.e. E1 + E2) by concatenating their data arrays 
%       along the first (temporal) dimension. E1 and E2 must have the same
%       number of channels.
%
%   [Eph,tIdx] = Eph.range(tLim) returns a new object whose data array 
%       contains only the samples that fall between time limits 
%       [tLim(1), tLim(2)]. Also returns the logical vector (tIdx) whose 
%       "true" indices correspond to the samples extracted from the 
%       original data array.
%
%   Eph.chan(idx) returns a new object whose data array contains only the
%       channel indices (idx) from the original object
%
%   [Eph,trialNo] = Eph.trial(num) returns a new object whose data array
%       contains only the trial numbers (num) from the original object. num
%       can also be specified as 'rand', in which case a random trial
%       number will be drawn and returned (trialNo)
%
%   [Eph,v] = Eph.cat_trials ...
%
%   Eph.rect rectifies (i.e. takes absolute value) the data array
%
%   Eph.filt(...) is a flexible method for filtering the data
%       filt(type,n,Wn) filters with an nth-order Butterworth filter with
%           passband Wn in Hz. If type is 'bandpass' and Wn a 2-element
%           vector, filt designs a bandpass filter. If type is 'low' or
%           'high' and Wn a scalar, filt designs a lowpass or highpass
%           filter, respectively.
%       filt('diff',n) applies an nth-order differentiating filter
%       filt('gau',sd) smooths with a Gaussian filter with standard
%           deviation sd in seconds
%       filt('box',wid) smooths with a boxcar filter of width wid in
%           seconds.
%
%   [Eph,sigma] = Eph.normalize returns a new object whose data array is
%       normalized to the estimated standard deviation of the noise on each
%       channel (sigma).
%
%   Eph.residual(Neu) ...
%
%   [Spk,idx] = Eph.detect_spikes detects spikes using FINDPULSES and
%       returns spike indices (idx) in a Spike object (Spk). Accepts all
%       optional and variable arguments as FINDPULSES
%
%   Eph.plot(...) plots all data. If array contains multiple channels,
%       plots as a set of vertically stacked subplots with linked x-axes.
%       plot(...,'fig',code) clears current axes if code is 'clf' or opens
%           a new figure if code is 'new'
%       plot(...,'sigma',sig) overlays a horizontal line at y = +/- sig
%           (useful for visualizing the threshold set by particular
%           multiples of the noise standard deviation)
%       plot(...,'Spk',S) overlays markers at the indices contained in the 
%           Spike object, S
%       plot(...,'Neu',N) overlays waveforms for each unit contained in the
%           Neuron object, N
%
%
% EXAMPLE(S) 
%
%
% IMPLEMENTATION
% Other m-files required: none
% Subfunctions: FINDPULSES
% MAT-files required: none
%
% SEE ALSO: FINDPULSES

% Authors: Najja Marshall
% Emails: njm2149@columbia.edu
% Dated:

classdef Ephys
    properties (Access = public)
        Fs = 1e3
        data
        alignIndex = 1
        chanNo
    end
    properties (SetAccess = private)
        time
        nSamples
        nChannels
        nTrials
        trialNo
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
        function [obj,trialNo] = trial(obj,num)
            % selects specified trials from data
            if isscalar(num)
                assert(num>=1 && num<=obj.nTrials)
                trialNo = num;
            elseif ischar(num) && strcmp(num,'rand')
                trialNo = randi(obj.nTrials);
            else
                error('Provide a scalar trial number or specify as ''rand''.')
            end
            obj.data = obj.data(:,:,trialNo);
            obj.trialNo = trialNo;
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
                        Wn = arg2;
                        
                        [b,a] = butter(n, Wn/(obj(ii).Fs/2), type);
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
%         function [res,vEst] = residual(obj,Neu)
%             % compute residual
%             assert(isa(Neu,'Neuron'))
%             vEst = zeros(obj.nSamples, obj.nChannels);
%             idxLim = [1 length(obj.time)];
%             for ch = 1:obj.nChannels
%                 for un = 1:Neu.nUnits
%                     si = find(Neu.spikes(:,un)); % + obj.alignIndex;
%                     si = si(si>idxLim(1) & si<idxLim(2));
%                     x = zeros(obj.nSamples,1);
%                     x(si) = 1;
%                     vEst(:,ch) = vEst(:,ch) + conv(x,Neu.waveform(:,ch,un),'same');
%                 end
%             end
%             eval(sprintf('yEst = %s(yEst);',class(obj.data)))
%             res = obj.data - vEst;
%         end
        % -----------------------------------------------------------------
        % WRAPPER FUNCTIONS
        % -----------------------------------------------------------------
        function [Spk,idx] = detect_spikes(obj,varargin)
            % detect spike events and return a Spike object
            isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;
            P = inputParser;
            addParameter(P, 'thresh', 4, @(x) isscalarnum(x,-Inf,Inf))
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
                Spk(ii).nSamples = obj.nSamples;
            end
        end
        % -----------------------------------------------------------------
        % PLOTTING
        % -----------------------------------------------------------------
        function plot(obj,varargin)
            % plot
            P = inputParser;
            addParameter(P,'fig',[],@(x) ischar(x) && ismember(x,{'clf','new'}))
            addParameter(P,'sigma',[],@isnumeric)
            addParameter(P,'Spk',[],@(x) isa(x,'Spike'))
            addParameter(P,'Neu',[],@(x) isa(x,'Neuron'))
            addParameter(P,'txt',[],@ischar)
            parse(P,varargin{:})
            
            obj.set_fig(P.Results.fig);
            
            ax = [];
            for ii = 1:obj.nChannels
                if obj.nChannels > 1
                    ax(ii) = subplot(obj.nChannels,1,ii);
                end
                plot(obj.time,obj.data(:,ii),'k')
                if ~isempty(P.Results.sigma)
                    hold on
                    plot(get(gca,'xlim'),P.Results.sigma*[1 1],'r')
                    plot(get(gca,'xlim'),-P.Results.sigma*[1 1],'r');
                end
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
        function fh = set_fig(~,code)
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