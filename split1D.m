%% SPLIT1D unsupervised clustering via 1-dimensional splits
% Function details
%
% SYNTAX
%   outputs = functiontemplate(inputs, varargin)
%
% REQUIRED INPUTS
%   reqIn (class): description
%
% OPTIONAL INPUTS
%   optIn (class): description
%
% PARAMETER INPUTS
%   'parameterName', <argument class>: description (default: )
%
% OUTPUTS
%   out1 (class): description
%
% EXAMPLE(S) 
%
% % 
% r = [.5+(rand(1e4,1)*2-.5); 8+(rand(1e4,1)*2-.5)];
% theta = [rand(1e4,1)*2*pi; rand(1e4,1)*pi];
% X = r.*[cos(theta),sin(theta)];
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

function label = split1D(X, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'SPLIT1D';

% validation functions
% isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'X', @isnumeric)
% addOptional(P, 'optIn', default, validationFunction)
addParameter(P, 'alpha', 1.2, @isscalar)
addParameter(P, 'maxSplit', Inf, @isscalar)
addParameter(P, 'maxRound', 20, @isscalar)
addParameter(P, 'pTest', 0.9, @isscalar)
addParameter(P, 'pThresh', 0.001, @isscalar)
addParameter(P, 'featureMaps', {'iso'}, @(x) iscell(x) && all(ismember(x,{'iso','sphere','rbf'})))

% clear workspace (parser object retains the data while staying small)
parse(P, X, varargin{:});
clear ans % varargin

%%

[nObs,nFeat] = size(X);
label = ones(nObs,1);

for iFeat = 1:length(P.Results.featureMaps)
    
    % unique labels
    uqLab = unique(label);
    trySplit = true(size(uqLab));
    noSplit = [];
    nSplits = 0;
    
    while any(trySplit) && nSplits<P.Results.maxSplit
        
        % working label
        labNo = uqLab(find(trySplit,1));
        lid = find(label==labNo);
        nl = length(lid);
        
        if nl <= nFeat
            trySplit(labNo) = false;
            noSplit = [noSplit, labNo];
            continue
        end
        
        Y = X(lid,:);
        
        switch P.Results.featureMaps{iFeat}
            
            case 'sphere'
                
                % mean-center
                Y = Y - mean(Y,1);
                
                % whiten
                [pcs,~,eigenvals] = pca(Y);
                Y = Y*pcs*diag(eigenvals)^(-1/2);
                
                % convert to spherical coordinates
                if nFeat == 2
                    Y = [sqrt(sum(Y.^2,2)), atan2(Y(:,2),Y(:,1))];
                    
                elseif nFeat == 3
                    Y = [sqrt(sum(Y.^2,2)), atan2(Y(:,2),Y(:,1)), atan2(sqrt(sum(Y(:,1:2).^2,2)),Y(:,3))];
                end
                
            case 'rbf'
                
                % mean-center
                Y = Y - mean(Y,1);
                
                % whiten
                [pcs,~,eigenvals] = pca(Y);
                Y = Y*pcs*diag(eigenvals)^(-1/2);
                
                % augment with RBF
                Y = [Y, sqrt(sum(Y.^2,2))];
        end
        
        % project data into 1D space
        [~,w,Ywht] = ica(Y,'contrast','tanh','maxRound',P.Results.maxRound);
        Ywht = Ywht';
        w = w';
        w = [round(w,4);eye(size(Ywht,2))];
        jj = 1;
        while jj <= size(w,1)
            w = setdiff(w,-w(jj,:),'rows');
            jj = jj+1;
        end
        w = unique(w,'rows')';
        
        % get spacings between ordered projections
        yProj = Ywht*w;
        yProjSrt = sort(yProj,1);
        dx = diff(yProjSrt,[],1);
        
        [pVal,cutPt] = deal(zeros(1,size(w,2)));
        for ii = 1:size(w,2)
            
            if isfinite(P.Results.pThresh)
                % Hartigan's dip test
                [~,pVal(ii)] = HartigansDipSignifTest(yProjSrt(:,ii), 1/P.Results.pThresh);
            else
                pVal(ii) = intmax;
            end
            
            % only fit central part of spacings
            minCut = ceil(nl*(1-P.Results.pTest)/2);
            [~,cutPt(ii)] = max(dx(minCut:end-minCut+1,ii));
        end
        
        pCrit = pVal<P.Results.pThresh;
        
        if any(pCrit)
            
            cutPt = cutPt+minCut-1;
            
            if nnz(pCrit) > 1
                [~,bestProj] = max(diag(dx(cutPt,:))' .* pCrit);
            else
                [~,bestProj] = min(pVal);
            end
            
            % new labels
            lNew = ones(nl,1);
            lNew(yProj(:,bestProj) > mean(yProjSrt([0 1]+cutPt(bestProj),bestProj))) = 2;
            
            if nnz(lNew==1)<=nFeat || nnz(lNew==2)<=nFeat
                trySplit(labNo) = false;
                noSplit = [noSplit, labNo];
            else
                label(lid) = 0*label(lid);
                newID = setdiff(1:2+max(label),unique(label));
                newID = newID(1:2);
                label(lid(lNew==1)) = newID(1);
                label(lid(lNew==2)) = newID(2);
                
                uqLab = unique(label);
                trySplit = true(size(uqLab));
                trySplit(ismember(uqLab,noSplit)) = false;
                
                nSplits = nSplits+1;
            end
        else
            trySplit(labNo) = false;
            noSplit = [noSplit, labNo];
        end
    end
end