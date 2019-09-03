%% WAVETEMPLATE waveform template
% Constructs the multi-channel waveform template from a set of waveforms
%
% SYNTAX
%   [u, srtIdx] = wavetemplate(w, varargin)
%
% REQUIRED INPUTS
%   w (cell): array of waveforms of dimensions [units x channels]
%
% OPTIONAL INPUTS: none
%
% PARAMETER INPUTS
%   'sort', <logical>: if true, sorts set of templates by their average
%       energy across channels (default: false)
%
%   'plot', <logical>: if true, plots templates (default: false)
%
%   'method', <string>: if 'mean' (default), estimates template as the mean
%       over all waveform shapes. If 'pmax', estimates template as the
%       smoothed maximum of the probability density at each sample point
%
%   'Fs', <scalar>: sample rate in Hz (default: 3e4)
%
%   'frame', <scalar>: fraction of waveform duration used to compute
%       templates (default: 1)
%
% OUTPUTS
%   u (numeric): 3D multi-channel template array of dimensions 
%       [wave length x channels x units]
%
%   srtIdx (numeric): vector of indices corresponding to the sort order (if
%       applied)
%
% EXAMPLE(S) 
%
%
% IMPLEMENTATION
% Other m-files required: PLOTWAVETEMPLATE
% Subfunctions: PLOTWAVETEMPLATE
% MAT-files required: none
%
% SEE ALSO: PLOTWAVETEMPLATE

% Authors: Najja Marshall
% Emails: njm2149@columbia.edu
% Dated: February 2019

function [u, srtIdx] = wavetemplate(w, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'WAVETEMPLATE';

% validation functions
isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'w', @iscell)
addParameter(P, 'sort', false, @islogical)
addParameter(P, 'plot', false, @islogical)
addParameter(P, 'method', 'mean', @(x) ischar(x) && ismember(x,{'mean','pmax'}))
addParameter(P, 'Fs', 3e4, @(x) isscalarnum(x,0,Inf))
addParameter(P, 'frame', 1, @(x) isscalarnum(x,0-eps,1+eps))

% clear workspace (parser object retains the data while staying small)
parse(P, w, varargin{:});
clear ans varargin

%%

[nUnits,nChannels] = size(w);
waveLen = size(w{1},2);

% set frame
len = round(waveLen*P.Results.frame);
win = waveLen/2 + (1-len/2:len/2);

% average and reshape to length x channel x unit
u = cellfun(@(x) mean(x(:,win),1)',w,'uni',false);

if strcmp(P.Results.method, 'mean')
    
    u = cell2mat(reshape(u',1,nChannels,nUnits));
    
else % re-estimate template based on maximum probability density
        
    uMax = cellfun(@max,u);
    uMin = cellfun(@min,u);
    u = cell2mat(reshape(u',1,nChannels,nUnits));
    
    for ii = 1:nUnits
        for jj = 1:nChannels
            yl = [min(-500,2*uMin(ii,jj)), max(500,2*uMax(ii,jj))];
            
            [~,b,lc] = histfun(w{ii,jj},'lim',yl);
            [~,maxIdx] = max(lc,[],1);
            u(:,jj,ii) = smooth1D(b(maxIdx)',P.Results.Fs,'gau','sd',1e-4);
        end 
    end
end

% sort by energy
if P.Results.sort
    ener = squeeze(sum(u.^2,1));
    ener = mean(ener,1);
    [~,srtIdx] = sort(ener);
    u = u(:,:,srtIdx);
else
    srtIdx = 1:nUnits;
end

% plot
if P.Results.plot
    plotwavetemplate(u);
end