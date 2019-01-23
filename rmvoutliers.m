%% RMVOUTLIERS remove outliers
% Function details
%
% SYNTAX
%   outputs = functiontemplate(inputs, varargin)
%
% REQUIRED INPUTS
%   reqIn (class) - description
%
% OPTIONAL INPUTS
%   optIn (class) - description
%
% PARAMETER INPUTS
%   'parameterName' <class> - description (default: )
%
% OUTPUTS
%   out1 (class) - description
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

function newLabel = rmvoutliers(X, label, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'RMVOUTLIERS';

% validation functions
isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;
isint = @(x,lb,ub) isscalar(x) && x==round(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'X', @isnumeric)
addRequired(P, 'label', @isnumeric)
addParameter(P, 'clus', [], @isnumeric)
addParameter(P, 'precision', 2, @(x) isint(x,0,5))
addParameter(P, 'cut', 0.95, @(x) isscalarnum(x,0,1))

% clear workspace (parser object retains the data while staying small)
parse(P, X, label, varargin{:});
clear ans varargin

%%

% input dimensions
T = size(X,2);

% unique clusters
unqClus = P.Results.clus;
if isempty(unqClus)
    unqClus = setdiff(unique(label),min(label):0);
end
nClus = length(unqClus);

% copy input labels
newLabel = label;

% loop clusters
for ii = 1:nClus
    
    clusIdx = find(label==unqClus(ii));
    
    if length(clusIdx) > 1
        
        y = X(clusIdx,:);
        
        % 2D histogram
        [counts,bins] = histfun(y);
        nBins = length(bins);
        counts = mat2cell(counts, nBins, ones(1,T));
        
        % filter counts
        filtCounts = cellfun(@(x) smooth1D(x,1,'gau','sd',25), counts, 'uni', false);
        
        % threshold
        normCounts = cell2mat(filtCounts)./cellfun(@sum,filtCounts);
        thr = max(normCounts,[],1);
        for jj = 1:T
            pow = floor(log10(thr(jj)));
            for kk = 0:P.Results.precision
                step = 10^(pow-kk);
                while sum(normCounts(:,jj).*(normCounts(:,jj)>thr(jj))) < P.Results.cut
                    thr(jj) = thr(jj)-step;
                end
                if kk < P.Results.precision
                    thr(jj) = thr(jj)+step;
                end
            end
        end
        
        % identify bounds above threshold
        normCounts = mat2cell(normCounts, nBins, ones(1,T));
        bounds = cellfun(@(x,thresh) bins([find(x>thresh,1), find(x>thresh,1,'last')]),...
            normCounts, num2cell(thr),'uni',false);
        
        % mark entries for removal that exceed bounds
        rmvIdx = false(size(y,1),1);
        for jj = 1:T
            rmvIdx = rmvIdx | y(:,jj)<bounds{jj}(1) | y(:,jj)>bounds{jj}(2);
        end
        
        % flip the sign on the label of outliers
        newLabel(clusIdx(rmvIdx)) = -label(clusIdx(rmvIdx));
    end
end