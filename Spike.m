%% Spike Data Class
%
% CONSTRUCTOR
%   Spk = Spike(Fs, index, waveform, label)
%
% PROPERTIES (public access)
%   Fs <scalar>: sampling frequency in Hz
%   index <numeric>: vector of spike indices
%   waveform <numeric>: array of waveform snippets extracted from the data
%       at each spike index
%   label <numeric>: vector of assignment labels for each spike event.
%       Unique, non-negative values correspond to distinct clusters.
%   nSamples <scalar>: number of data points in the full data array in
%       which spike events were detected
%   maxFeat <scalar>: maximum number of features stored in memory for each
%       waveform (default: 20)
%   clusDim <scalar>: number of features used for clustering (default: 3)
%
% PROPERTIES (private set access)
%   count: number of spike events
%   waveLength: number of samples in each waveform
%   feature: projections of the waveforms into the first maxFeat principal
%       components, inferred across all waveforms
%   template: waveform templates for each cluster
%   Clus: structure containing the number of dimensions used to perform the
%       clustering, the ID (label) of each cluster, size (number of
%       spikes), standard deviation at each time step, and entropy (for a
%       multivariate Gaussian)
%
% METHODS (public access)
%   Spk.indices equivalent to "index" if Spk is a scalar object. For an
%       object array, returns a cell array containing the spike indices for
%       each object. Further separates spike indices by cluster if "label"
%       contains any non-negative values.
%
%   Spk.sparsify converts spike indices into a sparse logical array of
%       dimensions nSamples x 1. If "label" contains any non-negative
%       values, will return an nSamples x n array where n is the number of
%       unique non-negative labels.
%
%   Spk.get_waveforms(v,waveDur) returns a new object with waveforms of
%       duration waveDur (seconds) extracted from the data v, centered at
%       each spike index. If Spk is an object array and v has the same
%       number of columns as the length of Spk, then the waveforms for each
%       object are obtained from each consecutive column in v.
%
%   Spk.align(v,[method],[waveDur],[refDur]) aligns spike waveforms using
%       user-specified method. If unspecified, method defaults to 'asym',
%       which uses an asymetric wavelet as the target template for 
%       alignment. If method is 'sym', uses a symmetric wavelet as the
%       target template. If method is 'self', aligns all waveforms within 
%       each cluster to the first waveform. If method is 'peak', uses a
%       peak-based alignment function (ALIGNWAVEFORMS). If using 'peak',
%       must also specify a waveform alignment duration (waveDur) and a 
%       refractory duration (refDur), both in units of seconds. Spike
%       indices that coincide within refDur seconds are removed.
%
%   Spk.cluster(nPC,maxClus) clusters waveforms using a Dirichlet-Process
%       Gaussian Mixture Model with nPC features and a truncation of
%       maxClus on the number of allowable latent clusters. Calls
%       MIXGAUSSVB to perform variational inference.
%
%
%
% EXAMPLE(S) 
%
%
% IMPLEMENTATION
% Other m-files required: none
% Subfunctions: ALIGNWAVEFORMS, RMVOUTLIERS, SMOOTH1D, STA, HISTFUN, PLOTMARKERS, MIXGAUSSVB
% MAT-files required: none
%
% SEE ALSO: ALIGNWAVEFORMS

% Authors: Najja Marshall
% Emails: njm2149@columbia.edu
% Dated:

classdef Spike
    properties
        Fs = 1e3
        index
        waveform
        label
        nSamples
        maxFeat = 20
        clusDim = 3
    end
    properties (SetAccess = private)
        count = 0
        waveLength
        feature
        template
        Clus
    end
    methods
        % -----------------------------------------------------------------
        % CONSTRUCTOR
        % -----------------------------------------------------------------
        function obj = Spike(Fs, index, waveform, label)
            if nargin > 0 && isnumeric(Fs) && isscalar(Fs)
                obj.Fs = Fs;
            end
            if nargin > 1 && isnumeric(index)
                obj.index = index(:);
            end
            if nargin > 2 && isnumeric(waveform)
                obj.waveform = waveform;
            end
            if nargin > 3 && isnumeric(label)
                obj.label = label(:);
            end
        end
        % -----------------------------------------------------------------
        % SET METHODS
        % -----------------------------------------------------------------
        function obj = set.index(obj,index)
            if isnumeric(index)
                obj.index = index;
                obj = obj.update_count;
            end
        end
        function obj = set.waveform(obj,wave)
            if isnumeric(wave)
                obj.waveform = wave;
                obj = obj.update_wavelen;
                obj = obj.update_template;
                obj = obj.update_feature;
            end
        end
        function obj = set.label(obj,label)
            if isnumeric(label)
                obj.label = label(:);
                obj = obj.update_template;
                obj = obj.update_clus;
            end
        end
        function obj = set.maxFeat(obj,maxFeat)
            if isnumeric(maxFeat) && isscalar(maxFeat)
                obj.maxFeat = maxFeat;
            end
        end
        % -----------------------------------------------------------------
        % GET METHODS
        % -----------------------------------------------------------------
        function lab = get.label(obj)
            if isempty(obj.label) && ~isempty(obj.waveform)
                lab = zeros(size(obj.waveform,1),1,'uint8');
            else
                lab = obj.label;
            end
        end
        function nPts = get.nSamples(obj)
            if isempty(obj.nSamples)
                if isempty(obj.index)
                    nPts = 0;
                else
                    nPts = max(obj.index);
                end
            else
                nPts = obj.nSamples;
            end
        end
        function feat = get.feature(obj)
            feat = obj.feature;
        end
        function nSpk = counts(obj)
            if length(obj) == 1
                nSpk = obj.count;
            else
                nSpk = arrayfun(@(x) x.count, obj);
            end
        end
        % -----------------------------------------------------------------
        % INDEX MANAGEMENT
        % -----------------------------------------------------------------
        function idx = indices(obj)
            % return all spike indices
            idx = cell(length(obj),1);
            for ii = 1:length(obj)
                if nnz(obj(ii).label)>0
                    unqLabels = unique(obj(ii).label(obj(ii).label>0));
                    idx{ii} = cellfun(@(n) obj(ii).index(obj(ii).label==n), num2cell(unqLabels), 'uni', false);
                else
                    idx{ii} = obj(ii).index;
                end
            end
            if length(obj) == 1
                idx = cat(1,idx{:});
            end
        end
        function s = sparsify(obj)
            % return spike locations as a sparse array
            assert(isequal(obj.nSamples,obj(1).nSamples),'All spike objects must have the same number of samples')
            s = cell(1,length(obj));
            for ii = 1:length(obj)
                unitNo = unique(obj(ii).label(obj(ii).label>0));
                nUnits = length(unitNo);
                if nUnits == 0
                    s{ii} = sparse(obj(ii).index,1,true,obj(ii).nSamples,1);
                else
                    s{ii} = spalloc(obj(ii).nSamples,nUnits,obj(ii).count);
                    for jj = 1:nUnits
                        idx = obj(ii).index(obj(ii).label==unitNo(jj));
                        s{ii}(:,jj) = sparse(idx,1,true,obj(ii).nSamples,1);
                    end
                end
            end
            s = logical(cell2mat(s));
        end
        % -----------------------------------------------------------------
        % WAVEFORM MANAGEMENT
        % -----------------------------------------------------------------
        function obj = get_waveforms(obj,v,waveDur)
            % get spike waveforms via spike-triggered average
            assert(size(v,2)==length(obj),'v must have same number of columns as entries in the object array')
            waveLen = round(waveDur*obj(1).Fs);
            waveLen = waveLen + mod(waveLen,2);
            for ii = 1:length(obj)
                [w,s] = spktrig(v(:,ii),num2cell(obj(ii).index),waveLen);
                obj(ii).waveform = cell2mat(w);
                obj(ii).index = cell2mat(s);
            end
        end
        function obj = align(obj,v,method,waveDur,refDur)
            
            if nargin == 2
                method = 'asym';
                waveDur = [];
                refDur = [];
            end
            
            switch method
                
                case 'self' % align within cluster
                    
                    for ii = 1:length(obj)
                        
                        len = 2*obj(ii).waveLength-1;
                        nPad = 2^nextpow2(len);
                        
                        nClus = length(obj(ii).Clus.id);
                        for jj = 1:nClus
                            
                            cid = obj(ii).label == obj(ii).Clus.id(jj);
                            
                            w = obj(ii).waveform(cid,:)';
                            
                            xc = ifft(fft(flipud(w(:,1)),nPad) .* fft(w,nPad,1),[],1);
                            [~,pkIdx] = max(xc);
                            obj(ii).index(cid) = obj(ii).index(cid) - (pkIdx(:)-obj(ii).waveLength);
                            
                            [w,s] = spktrig(v(:,ii),num2cell(obj(ii).index(cid)),obj(ii).waveLength);
                            obj(ii).waveform(cid,:) = cell2mat(w);
                            obj(ii).index(cid) = cell2mat(s);
                        end
                    end
                    
                case 'peak' % peak-based
                    
                    for ii = 1:length(obj)
                        [obj(ii).waveform,obj(ii).index] = alignwaveforms(v(:,ii),obj(ii).Fs,obj(ii).index,'waveDur',waveDur,'refDur',refDur);
                    end
                    
                otherwise
                    
                    waveLen = max(arrayfun(@(x) x.waveLength,obj));
                    len = 2*waveLen-1;
                    nPad = 2^nextpow2(len);
                    
                    if strcmp(method,'sym') % symmetric wavelet
                        
                        filtBank = cwtfilterbank('SignalLength',waveLen,...
                            'SamplingFrequency',obj(1).Fs,...
                            'waveletParameters',[1 3],...
                            'FrequencyLimits',fliplr(round(1e3./[.1 10])));
                        wvMorse = wavelets(filtBank);
                        wvTmp = real(wvMorse(30,:));
                        
                    else % asymmetric wavelet
                        
                        waveParams = [1 3; 3 9];
                        for ii = 1:2
                            filtBank = cwtfilterbank('SignalLength',waveLen,...
                                'SamplingFrequency',obj(1).Fs,...
                                'waveletParameters',waveParams(ii,:),...
                                'FrequencyLimits',fliplr(round(1e3./[.1 10])));
                            wvMorse = wavelets(filtBank);
                            y = real(wvMorse(round(size(wvMorse,1)/2),:));
                            
                            if ii == 1
                                wvTmp = y;
                            else
                                stitchPt = 1+waveLen/2;
                                wvTmp(stitchPt:end) = (wvTmp(stitchPt)/y(stitchPt)) * y(stitchPt:end);
                            end
                        end
                    end
                    wvTmp = fft(wvTmp',nPad,1);
                    
                    spkLim = [1+waveLen/2, obj(1).nSamples-waveLen/2+1];
                    for ii = 1:length(obj)
                        if isempty(obj(ii).index)
                            continue
                        end
                        
                        pAligned = 0;
                        dp = Inf;
                        while pAligned < 0.99 && dp > 1e-3
                            xc = ifft(wvTmp .* fft(flipud(obj(ii).waveform'),nPad,1),[],1);
                            [~,pkIdx] = max(xc);
                            oldIdx = obj(ii).index;
                            obj(ii).index = obj(ii).index - (pkIdx(:)-obj(ii).waveLength);
                            
                            dp = nnz(oldIdx==obj(ii).index)/length(obj(ii).index) - pAligned;
                            pAligned = nnz(oldIdx==obj(ii).index)/length(obj(ii).index);
                            
                            outOfBounds = obj(ii).index < spkLim(1) | obj(ii).index > spkLim(2);
                            obj(ii).index(outOfBounds) = [];
                            obj(ii).label(outOfBounds) = [];
                            obj(ii).waveform(outOfBounds,:) = [];
                            
                            [w,s] = spktrig(v(:,ii),num2cell(obj(ii).index),obj(ii).waveLength);
                            obj(ii).waveform = cell2mat(w);
                            obj(ii).index = cell2mat(s);
                        end
                        
                        % remove refractory period violations
                        refViol = diff(obj(ii).index/obj(ii).Fs) < 5e-4;
                        obj(ii).index(refViol) = [];
                        obj(ii).label(refViol) = [];
                        obj(ii).waveform(refViol,:) = [];
                        
                        % center mean waveform mass
                        wEnv = smooth1D(abs(obj(ii).waveform),obj(ii).Fs,'gau','sd',5e-4,'dim',2);
                        [~,pkIdx] = max(mean(wEnv,1));
                        ds = obj(ii).waveLength/2 - pkIdx;
                        [w,s] = spktrig(v(:,ii),num2cell(obj(ii).index - ds),obj(ii).waveLength);
                        obj(ii).waveform = cell2mat(w);
                        obj(ii).index = cell2mat(s);
                    end
            end
        end
        % -----------------------------------------------------------------
        % CLUSTER
        % -----------------------------------------------------------------
        function obj = cluster(obj,nPC,maxClus)
            % cluster waveforms via DP-GMM
            for ii = 1:length(obj)
                if obj(ii).count == 0
                    continue
                end
                obj(ii).clusDim = nPC;
                obj(ii).label = mixGaussVb(obj(ii).feature(:,1:nPC)',maxClus)';
            end
        end
        % -----------------------------------------------------------------
        % LABEL MANAGEMENT
        % -----------------------------------------------------------------
        % set
        function obj = set_labels(obj,allLabels,grp)
            lNew = zeros(obj.count,1,'uint8');
            for ii = 1:length(grp)
                lNew(ismember(allLabels,grp{ii})) = ii;
            end
            obj.label(lNew>0) = lNew(lNew>0) + max(obj.label);
        end
        % remove select labels
        function obj = rmv_labels(obj,labelNo)
            obj.label(ismember(abs(obj.label),labelNo)) = 0;
        end
        % clear all labels
        function obj = clr_labels(obj,varargin)
            P = inputParser;
            addParameter(P,'except',[],@isnumeric)
            parse(P,varargin{:})
            for ii = 1:length(obj)
                rmvNo = setdiff(unique(abs(obj(ii).label)),P.Results.except);
                obj(ii).label(ismember(abs(obj(ii).label),rmvNo)) = 0;
            end
        end
        % reset numbering
        function obj = reset_labels(obj)
            newLabels = zeros(obj.count,1,'uint8');
            unqLabels = unique(obj.label);
            unqLabels = unqLabels(unqLabels>0);
            for ii = 1:length(unqLabels)
                newLabels(obj.label ==  unqLabels(ii)) =  ii;
                newLabels(obj.label == -unqLabels(ii)) = -ii;
            end
            obj.label = newLabels;
        end      
        % triage (remove outliers)
        function obj = triage(obj,cutoff)
            if length(cutoff) == 1
                cutoff = cutoff*ones(1,length(obj));
            end
            for ii = 1:length(obj)
                if ~isempty(obj(ii).label)
                    obj(ii).label = rmvoutliers(double(obj(ii).waveform),obj(ii).label,'cut',cutoff(ii));
                end
            end
        end
        % merge
        function obj = mrg_labels(obj,labelGrp)
            for ii = 1:length(labelGrp)
                % shift spike indices
                unitLabel = setdiff(obj.label,min(obj.label):0);
                wT = obj.template(:,unitLabel==labelGrp{ii}(1),1);
                for jj = 2:length(labelGrp{ii})
                    lid = obj.label==labelGrp{ii}(jj);
                    spkIdx = obj.index(lid);
                    w = double(obj.waveform(lid,:))';
                    for kk = 1:length(spkIdx)
                        [c,l] = xcorr(wT,w(:,kk));    
                        [~,maxIdx] = max(c);
                        spkIdx(kk) = spkIdx(kk)+l(maxIdx);
                    end
                    obj.index(lid) = spkIdx;
                end
                obj.label(ismember(obj.label,labelGrp{ii})) = labelGrp{ii}(1);
            end
            newLabels = zeros(obj.count,1,'uint8');
            unqLabels = unique(obj.label);
            unqLabels = unqLabels(unqLabels>0);
            for ii = 1:length(unqLabels)
                newLabels(obj.label==unqLabels(ii)) = ii;
            end
            obj.label = newLabels;                
        end
        % assign unsorted
        function obj = assign_unsorted(obj,thresh)
            if length(thresh) == 1
                thresh = thresh*ones(1,length(obj));
            end
            for ii = 1:length(obj)
                if isempty(obj(ii).label)
                    continue
                end
                % unsorted index
                unsrtIdx = find(obj(ii).label==0);
                nUnsrt = length(unsrtIdx);
                % unsorted waveforms
                waveUnsrt = double(obj(ii).waveform(unsrtIdx,:)');
                % compute similarity w.r.t. all templates
                nTemplates = size(obj(ii).template,2);
                wSim = zeros(nUnsrt,nTemplates);
                for jj = 1:nTemplates
                    % update similarity function
                    wTmp = obj(ii).template(:,jj,1);
                    tempVar = var(wTmp);
                    fn = @(x) 1 - ((x-wTmp)'*(x-wTmp))/(obj(ii).waveLength*tempVar);
                    for kk = 1:nUnsrt
                        wSim(kk,jj) = fn(waveUnsrt(:,kk));
                    end
                end
                [wSim,nearestTemplate] = max(wSim,[],2);
                obj(ii).label(unsrtIdx(wSim>thresh(ii))) = obj(ii).Clus.id(nearestTemplate(wSim>thresh(ii)));
            end
        end
        function chanNo = unit_channel(obj)
            % return channel numbers for each cluster on each channel
            allChanNo = num2cell(1:length(obj))';
            nUnits = arrayfun(@(S) size(S.template,2),obj,'uni',false);
            chanNo = cellfun(@(cNo,nU) repmat(cNo,nU,1), allChanNo, nUnits,'uni',false); 
            chanNo = cat(1,chanNo{:});
        end
        % -----------------------------------------------------------------
        % PLOTTING
        % -----------------------------------------------------------------
        
        % 2D histogram
        function plot_hist(obj,varargin)
            P = inputParser;
            addParameter(P,'mask',[],@islogical)
            addParameter(P,'label',[],@isnumeric)
            addParameter(P,'lim',[],@(x) length(x)==2 && x(2)>x(1))
            addParameter(P,'nBins',[],@(x) isnumeric(x))
            addParameter(P,'ppp',[],@isnumeric)
            addParameter(P,'sort','energy',@(x) ischar(x) && ismember(x,{'id','size','energy','entropy'}))
            parse(P,varargin{:})
            
            if length(obj) > 1
                figure
                nCol = ceil(sqrt(length(obj)));
                nRow = ceil(length(obj)/nCol);
                for ii = 1:length(obj)
                    if obj(ii).count > 0
                        subplot(nRow,nCol,ii)
                        histfun(obj(ii).waveform,'plot',true,'lim',P.Results.lim,'nBins',P.Results.nBins);
                    end
                end
            else
                
                % mask
                mask = P.Results.mask;
                if isempty(mask)
                    mask = true(size(obj.waveform,1),1);
                end
                
                if nnz(unique(obj.label)>0) < 2
                    mask = (obj.label == max(unique(obj.label)));
                    histfun(obj.waveform(mask,:),'plot',true,'lim',P.Results.lim,'nBins',P.Results.nBins);
                    title(sprintf('%i waveforms',nnz(mask)))
                    
                else
                    if isempty(P.Results.label)
                        unitNo = obj.Clus.id;
                    else
                        unitNo = P.Results.label;
                    end
                    nUnits = length(unitNo);
                    
                    [~,srtIdx] = sort(obj.Clus.(P.Results.sort));
                    
                    plotsPerPage = P.Results.ppp;
                    if isempty(plotsPerPage)
                        plotsPerPage = nUnits;
                    end
                    nPages = ceil(nUnits/plotsPerPage);
                    
                    nCol = ceil(sqrt(plotsPerPage));
                    nRow = ceil(plotsPerPage/nCol);
                    
                    for ff = 1:nPages
                        figure
                        for idx = 1:plotsPerPage
                            ii = idx+plotsPerPage*(ff-1);
                            if ii > nUnits
                                break
                            end
                            if nUnits > 1
                                subplot(nRow,nCol,idx)
                            end
                            if nnz(obj.label==unitNo(srtIdx(ii))) <= 1
                                continue
                            end
                            histfun(obj.waveform(obj.label==unitNo(srtIdx(ii)),:),'plot',true,'lim',P.Results.lim,'nBins',P.Results.nBins);
                            cid = obj.Clus.id == obj.Clus.id(srtIdx(ii));
                            title(sprintf('unit %i\nn=%i, H=%.2f, |w|^2=%.2f',unitNo(srtIdx(ii)),...
                                obj.Clus.size(cid), obj.Clus.entropy(cid), obj.Clus.energy(cid)))
                        end
                    end
                end
            end
        end
        
        % plot features
        function plot_feat(obj,featNo,colorLabels)
            if nargin == 2
                colorLabels = true;
            end
            
            % assemble feature values (0 corresponds to timestamp)
            nFeat = length(featNo);
            assert(nFeat==2 || nFeat==3, 'Specify 2 or 3 features')
            [featVal,axLabel] = deal(cell(nFeat,1));
            allFeatures = obj.feature;
            for ii = 1:nFeat
                if featNo(ii) == 0
                    featVal{ii} = obj.index/obj.Fs;
                    axLabel{ii} = 'timestamp (s)';
                else
                    featVal{ii} = allFeatures(:,featNo(ii));
                    axLabel{ii} = sprintf('PC %i',featNo(ii));
                end
            end
            
            if colorLabels
                % generate set of point markers
                colors = brewermap(8,'Set1');
                markers = {'.','*','o','+','x','s','d','^','v','>','<','p','h'};
                pm = reshape(repmat(markers, size(colors,1), 1),...
                    size(colors,1)*length(markers), 1);
                pm = [pm, repmat({'color'},size(pm,1),1), repmat(mat2cell(colors,ones(size(colors,1),1),3),length(markers),1)];
                % plot
                hold on
                clusNo = setdiff(unique(abs(obj.label)), min(obj.label):-1:0);
                nClus = length(clusNo);
                legText = [];
                % unlabeled
                if any(obj.label==0)
                    lid = obj.label==0;
                    if nFeat == 2
                        plot(featVal{1}(lid),featVal{2}(lid),'k.','color',0.8*ones(1,3))
                    else
                        plot3(featVal{1}(lid),featVal{2}(lid),featVal{3}(lid),'k.','color',0.8*ones(1,3))
                    end
                    legText = {'unsorted'};
                end
                % labeled
                if any(clusNo > 0)
                    for ii = 1:nClus
                        lid = obj.label==clusNo(ii);
                        if nFeat == 2
                            plot(featVal{1}(lid),featVal{2}(lid),pm{ii,:})
                        else
                            plot3(featVal{1}(lid),featVal{2}(lid),featVal{3}(lid),pm{ii,:})
                        end
                    end
                    legText = [legText, cellfun(@(n) sprintf('MU %i',n), num2cell(obj.Clus.id), 'uni', false)];
                end
                legend(legText)
            else
                if nFeat == 2
                    plot(featVal{1},featVal{2},'k.')
                else
                    plot3(featVal{1},featVal{2},featVal{3},'k.')
                end
            end
            
            % label axes
            xlabel(axLabel{1})
            ylabel(axLabel{2})
            if nFeat == 3
                zlabel(axLabel{3})
                view(3)
            end
        end 
        
        % templates
        function plot_templates(obj,varargin)
            P = inputParser;
            addParameter(P,'unit',[],@isnumeric)
            parse(P,varargin{:})
            
            lenObj = length(obj);
            
            nCol = ceil(sqrt(lenObj));
            nRow = ceil(lenObj/nCol);
            
            figure
            for ii = 1:lenObj
                if lenObj > 1
                    subplot(nRow,nCol,ii)
                end
                hold on
                
                unitNo = P.Results.unit;
                if isempty(unitNo)
                    unitNo = setdiff(unique(obj(ii).label),min(obj(ii).label):0);
                end
                nUnits = length(unitNo);
                if nUnits == 0
                    continue
                end
                
                cmap = brewermap(nUnits,'Spectral');
                
                unitCount = cellfun(@(x) nnz(obj(ii).label==x), num2cell(unitNo), 'uni', false);
                
                t = 1e3*(1:obj(ii).waveLength)/obj(ii).Fs;
                plot(t([1 end]),[0 0],'--','color',0.85*[1 1 1],'linewidth',0.5)
                fh = zeros(1,nUnits);
                for jj = 1:nUnits
                    fh(jj) = plot(t, obj(ii).template(:,jj,1), 'color', cmap(jj,:));
                    patch([t,fliplr(t)], [permute(sum(obj(ii).template(:,jj,:),3),[2 1 3]),...
                        fliplr(permute(-diff(obj(ii).template(:,jj,:),[],3),[2 1 3]))], cmap(jj,:),...
                        'EdgeAlpha',0, 'FaceAlpha',0.125)
                end
                legend(fh, cellfun(@(n,nc) sprintf('unit %i |%i|',n,nc), num2cell(unitNo), unitCount, 'uni', false))
                xlabel('time (ms)')
            end
        end
    end
    methods (Access = private)
        function obj = update_count(obj)
            obj.count = length(obj.index);
        end
        function obj = update_wavelen(obj)
            obj.waveLength = size(obj.waveform,2);
        end
        function obj = update_template(obj)
            clusNo = setdiff(unique(obj.label), min(obj.label):0);
            nClus = length(clusNo);
            if nClus > 0
                obj.template = zeros(obj.waveLength, nClus, 2);
                for ii = 1:nClus
                    w = obj.waveform(obj.label==clusNo(ii),:);
                    obj.template(:,ii,1) = mean(double(w),1);
                    obj.template(:,ii,2) = std(double(w),[],1)/sqrt(size(w,1)-1);                    
                end
            else
                obj.template = [];
            end
        end
        function obj = update_feature(obj)
            [~,obj.feature] = pca(double(obj.waveform));
            if size(obj.feature,2) > 0
                obj.maxFeat = min(obj.maxFeat, size(obj.feature,2));
                obj.feature = obj.feature(:,1:obj.maxFeat);
            end
        end
        function obj = update_clus(obj)
            if nnz(obj.label) == 0
                return
            end
            obj.Clus = struct('dim',{[]},'id',{[]},'size',{[]},'energy',{[]},'std',{[]},'entropy',{[]});
            obj.Clus.dim = obj.clusDim;
            obj.Clus.id = setdiff(unique(obj.label), min(obj.label):0);
            obj.Clus.id = obj.Clus.id(:)';
            [obj.Clus.size, obj.Clus.energy, obj.Clus.entropy] = deal(zeros(size(obj.Clus.id)));
            obj.Clus.std = zeros(length(obj.Clus.id),obj.waveLength);
            for ii = 1:length(obj.Clus.id)
                cid = obj.label == obj.Clus.id(ii);
                obj.Clus.size(ii) = nnz(cid);
                % waveforms
                obj.Clus.energy(ii) = sum(mean(double(obj.waveform(cid,:)),1).^2)/obj.waveLength;
                obj.Clus.std(ii,:) = std(double(obj.waveform(cid,:)),[],1);
                % features
                x = obj.feature(cid,1:obj.clusDim);
                mu = mean(x,1);
                S = (x-mu)'*(x-mu)/(nnz(cid)-1);
                obj.Clus.entropy(ii) = (obj.clusDim/2 + obj.clusDim/2*log(2*pi) + 1/2*log(det(S)))/log(2);
                if ~isreal(obj.Clus.entropy(ii)) || obj.Clus.entropy(ii) < 0
                    obj.Clus.entropy(ii) = Inf;
                end
            end
        end
    end
end
                
                
               
    
        