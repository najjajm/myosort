%% ICA independent components analysis
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

function [S,B,Z] = ica(X, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'ICA';

% validation functions
% isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'X', @isnumeric)
% addOptional(P, 'optIn', default, validationFunction)
addParameter(P, 'maxRound', 20, @isscalar)
addParameter(P, 'maxIter', 1e3, @isscalar)
addParameter(P, 'contrast', 'kurtosis', @(x) ischar(x) && ismember(x,{'kurtosis','tanh','exp'}))

% clear workspace (parser object retains the data while staying small)
parse(P, X, varargin{:});
clear ans varargin

%% Pre-processing

% ensure double
X = double(X);

% input dimensions
[~,nFeat] = size(X);

% mean center
muX = mean(X,1);
X = X-muX;

% whiten
[pcs,~,eigenvals] = pca(X);
Z = X*pcs*diag(eigenvals)^(-1/2);

% reshape
Z = Z';

% contrast functions
switch P.Results.contrast
    case 'kurtosis' % G = @(xx) xx.^4/4;
        g = @(xx) xx.^3;
        dg = @(xx) 3*xx.^2;   
        
    case 'tanh' % G = @(xx) log(cosh(xx));
        g = @(xx) tanh(xx);
        dg = @(xx) sech(xx).^2;
        
    case 'exp' % G = @(xx) -exp(-xx.^2/2);
        g = @(xx) xx.*exp(-xx.^2/2);
        dg = @(xx) exp(-xx.^2/2) .* (1-xx.^2); 
end

%% Projection pursuit

% initialize

B = [];
S = [];
contrast = zeros(P.Results.maxRound,P.Results.maxIter);

wStore = zeros(nFeat,P.Results.maxRound);

for ii = 1:P.Results.maxRound
    
    w = zeros(nFeat,P.Results.maxIter);
    w(:,1) = rand(nFeat,1);
    w(:,1) = w(:,1)/sqrt(w(:,1)'*w(:,1));

    for jj = 2:P.Results.maxIter
        
        % fixed point algorithm
        w(:,jj) = mean(Z.*g(w(:,jj-1)'*Z),2) - mean(dg(w(:,jj-1)'*Z))*w(:,jj-1);
        
        % orthonormalize
        if ~isempty(B)
            w(:,jj) = w(:,jj) - Borth*w(:,jj);
        end
        w(:,jj) = w(:,jj)/sqrt(w(:,jj)'*w(:,jj));
        
        % check convergence
        contrast(ii,jj) = abs(w(:,jj)'*w(:,jj-1) - 1);
        if abs(diff(contrast(ii,jj-1:jj))) < 1e-5
            break
        end
    end
    
    if ~isfinite(w(:,jj))
        break
    end
    
    % update orthogonalization matrix
    B = [B,w(:,jj)];
    Borth = B*B';
    
    % recover source estimate
    S = [S,Z'*w(:,jj)];
    
    wStore(:,ii) = w(:,jj);
end

