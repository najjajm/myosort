%% VIZ visualize sort progress

function viz(X,opt,sortBlock)
%% Load data

load([opt.savePath 'spikes'],'Spk')
load([opt.savePath 'labels'],'Lab')
load([opt.savePath 'templates'],'W')

%%

nChan = size(X,1);
waveLen = round(opt.Fs*opt.waveformDuration(2));

switch sortBlock
    
    case 'align'
        
        figure('Name','alignment')
        
        % plot initial detected waveforms
        subplot(121)
        hold on
        
        waveforms = myo.spktrig(X, double(Spk.detect), waveLen);
        nObs = size(waveforms,3);
        if nObs > 1000
            histfun(cellfun(@(x) squeeze(x)',mat2cell(waveforms,...
                ones(size(waveforms,1),1),size(waveforms,2),nObs),'uni',false),...
                'plot',true,'nBins',5e3);
        else
            myo.plotwaves(waveforms)
        end
        ylabel('au')
        title(sprintf('Detected\nn = %i',length(Spk.detect)))
        
        % plot aligned waveforms
        subplot(122)
        hold on
        
        waveforms = myo.spktrig(X, double(Spk.align), waveLen);
        nObs = size(waveforms,3);
        if nObs > 1000
            histfun(cellfun(@(x) squeeze(x)',mat2cell(waveforms,...
                ones(size(waveforms,1),1),size(waveforms,2),nObs),'uni',false),...
                'plot',true,'nBins',5e3);
        else
            myo.plotwaves(waveforms)
        end
        ylabel('au')
        title(sprintf('Aligned\nn = %i',length(Spk.align)))
        
    case 'cluster'
        
        waveforms = myo.spktrig(X, double(Spk.align), waveLen);
        nObs = size(waveforms,3);
        
        % plot waveforms in each cluster
        uql = unique(Lab.cluster);
        nClus = length(uql);
        unitIdx = [1:opt.maxSubplotsPerPage:nClus, 1+nClus];
        nFig = length(unitIdx)-1;
        
        for iFi = 1:nFig
            
            uNo = unitIdx(iFi):min(nClus,unitIdx(iFi+1)-1);
            if length(uNo)>1
                figure('Name',sprintf('clusters (units %i-%i)',uNo(1),uNo(end)))
            else
                figure('Name',sprintf('clusters (unit %i)',uNo))
            end
            
            for iUn = 1:length(uNo)
                
                subplot(1,opt.maxSubplotsPerPage,iUn)
                hold on
                
                lid = Lab.cluster==uNo(iUn);
                
                if nObs > 1000
                    histfun(cellfun(@(x) squeeze(x)',mat2cell(waveforms(:,:,lid),...
                        ones(size(waveforms,1),1),size(waveforms,2),nnz(lid)),'uni',false),...
                        'plot',true,'nBins',5e3);
                else
                    myo.plotwaves(waveforms(:,:,lid));
                end
                ylabel('au')
                title(sprintf('MU %i\nn = %i',uNo(iUn),nnz(lid)))
            end
        end
        
        % plot templates
        plotwavetemplate(W.cluster)
        set(gcf,'Name','templates')
        
    case 'merge'
        
        waveforms = myo.spktrig(X, double(Spk.merge), waveLen);
        nObs = size(waveforms,3);
        
        % plot waveforms in each cluster
        uql = unique(Lab.merge);
        nClus = length(uql);
        unitIdx = [1:opt.maxSubplotsPerPage:nClus, 1+nClus];
        nFig = length(unitIdx)-1;
        
        for iFi = 1:nFig
            
            uNo = unitIdx(iFi):min(nClus,unitIdx(iFi+1)-1);
            if length(uNo)>1
                figure('Name',sprintf('merged clusters (units %i-%i)',uNo(1),uNo(end)))
            else
                figure('Name',sprintf('merged clusters (unit %i)',uNo))
            end
            
            for iUn = 1:length(uNo)
                
                subplot(1,opt.maxSubplotsPerPage,iUn)
                hold on
                
                lid = Lab.merge==uNo(iUn);
                
                if nObs > 1000
                    histfun(cellfun(@(x) squeeze(x)',mat2cell(waveforms(:,:,lid),...
                        ones(size(waveforms,1),1),size(waveforms,2),nnz(lid)),'uni',false),...
                        'plot',true,'nBins',5e3);
                else
                    myo.plotwaves(waveforms(:,:,lid));
                end
                ylabel('au')
                title(sprintf('MU %i\nn = %i',uNo(iUn),nnz(lid)))
            end
        end
        
        % plot templates
        plotwavetemplate(W.merge)
        set(gcf,'Name','merged templates')
        
    case 'triage'
        
        waveforms = myo.spktrig(X, double(Spk.triage), waveLen);
        nObs = size(waveforms,3);
        
        % plot waveforms in each cluster
        uql = unique(Lab.merge);
        nClus = length(uql);
        unitIdx = [1:opt.maxSubplotsPerPage:nClus, 1+nClus];
        nFig = length(unitIdx)-1;
        
        for iFi = 1:nFig
            
            uNo = unitIdx(iFi):min(nClus,unitIdx(iFi+1)-1);
            if length(uNo)>1
                figure('Name',sprintf('triaged clusters (units %i-%i)',uNo(1),uNo(end)))
            else
                figure('Name',sprintf('triaged clusters (unit %i)',uNo))
            end
            
            for iUn = 1:length(uNo)
                
                subplot(1,opt.maxSubplotsPerPage,iUn)
                hold on
                
                lid = Lab.triage==uNo(iUn);
                
                if nObs > 1000
                    histfun(cellfun(@(x) squeeze(x)',mat2cell(waveforms(:,:,lid),...
                        ones(size(waveforms,1),1),size(waveforms,2),nnz(lid)),'uni',false),...
                        'plot',true,'nBins',5e3);
                else
                    myo.plotwaves(waveforms(:,:,lid));
                end
                ylabel('au')
                title(sprintf('MU %i\nn = %i',uNo(iUn),nnz(lid)))
            end
        end
        
        % plot templates
        plotwavetemplate(W.merge)
        set(gcf,'Name','triaged templates')
        
    case 'vet'
        
        plotwavetemplate(W.vet)
        set(gcf,'Name','vetted templates')
end

drawnow