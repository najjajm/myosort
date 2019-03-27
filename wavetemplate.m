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
% OUTPUTS
%   u (numeric): 3D multi-channel template array of dimensions 
%       [wave length x channels x units]
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

% add required, optional, and parameter-value pair arguments
addRequired(P, 'w', @iscell)
addParameter(P, 'sort', false, @islogical)
addParameter(P, 'plot', false, @islogical)

% clear workspace (parser object retains the data while staying small)
parse(P, w, varargin{:});
clear ans varargin

%%

[nUnits,nChannels] = size(w);

% average and reshape to length x channel x unit
u = cellfun(@(x) mean(x,1)',w,'uni',false);
u = cell2mat(reshape(u',1,nChannels,nUnits));

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
    plotwavetemplate(u)
end