%% MUTUALINFO mutual information
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

function [MI, normMI] = mutualinfo(spikes, Fs, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'MUTUALINFO';

% validation functions
isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'spikes', @(x) islogical(x) || iscell(x))
addRequired(P, 'Fs', @(x) isscalarnum(x,0,Inf))
addOptional(P, 'nT', [], @(x) isscalarnum(x,0,Inf))
addParameter(P, 'maxLag', 5e-3, @(x) isscalarnum(x,0,1))
addParameter(P, 'group', [], @(x) isnumeric(x) && (length(x)==size(spikes,2) || length(x)==length(spikes)))

% clear workspace (parser object retains the data while staying small)
parse(P, spikes, Fs, varargin{:});
clear ans varargin

%% 

grp = P.Results.group;
if islogical(spikes)
    % input dimensions
    [nDataPoints, nUnits] = size(spikes);
    
    % extract spike indices
    spikeIdx = cell(nUnits,1);
    for ii = 1:nUnits
        spikeIdx{ii} = find(spikes(:,ii));
    end
    
    if isempty(grp)
        grp = (1:size(spikes,2))';
    end
else
    nDataPoints = P.Results.nT;
    spikeIdx = spikes(:);
    if isempty(nDataPoints)
        warning('Data length unspecified. Using maximum spike index.')
        nDataPoints = max(cellfun(@max,spikes));
    end
    nUnits = length(spikeIdx);
    if isempty(grp)
        grp = (1:length(spikes))';
    end
end

% Shannon entropy
H = @(prob) -sum(arrayfun(@(p) min(0,p*log2(p)), prob));

% marginal (neural) entropies
nSpks = cellfun(@length, spikeIdx);
pNeu = nSpks/nDataPoints;
pNeu = [pNeu, 1-pNeu];
Hneu = cellfun(@(p) H(p), mat2cell(pNeu,ones(nUnits,1),2));

% time lags
maxLagSamp = round(Fs * P.Results.maxLag);
lags = num2cell(-maxLagSamp:maxLagSamp);
nLags = length(lags);

%% Mutual information

% lagged mutual information
MI = zeros(nUnits,nUnits,nLags);
for ii = 1:nUnits
    t0 = tic;
    for jj = 1:nUnits
        if jj>ii && grp(ii)~=grp(jj)
            
            % lagged spike indices
            sjLagged = cellfun(@(l) spikeIdx{jj}+l,lags,'uni',false);
            
            % joint spike probabilities
            nij = cellfun(@(sj) [...
                nnz( ismember(spikeIdx{ii},sj)),...
                nnz(~ismember(spikeIdx{ii},sj)),...
                nnz(~ismember(sj,spikeIdx{ii}))],sjLagged,'uni',false);
            pij = cellfun(@(x) x/nDataPoints, nij,'uni',false);
            pij = cellfun(@(p) [p, 1-sum(p)], pij,'uni',false);
            
            % joint entropy
            Hij = cellfun(@(p) H(p), pij);
            
            % mutual information
            MI(ii,jj,:) = Hneu(ii) + Hneu(jj) - reshape(Hij,1,1,nLags);
        end
    end
end

% normalized maximum mutual information
pkMI = max(MI,[],3);
normMI = zeros(size(pkMI));
for ii = 1:nUnits
    for jj = 1:nUnits
        if jj>ii && grp(ii)~=grp(jj)
            normMI(ii,jj) = 2*pkMI(ii,jj)/(Hneu(ii)+Hneu(jj));
        end
    end
end

% format outputs
if nUnits == 2
    MI = squeeze(MI(1,2,:));
    normMI = normMI(1,2);
else
    for kk = 1:nLags
        MI(:,:,kk) = MI(:,:,kk) + MI(:,:,kk)';
    end
    normMI = normMI + normMI';
end