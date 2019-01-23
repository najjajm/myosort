%% WAVETEMPLATE full function name
% Function details
%
% SYNTAX
%   outputs = functiontemplate(inputs, varargin)
%
% REQUIRED INPUTS
%   reqIn <class>: description
%
% OPTIONAL INPUTS
%   optIn <class>: description
%
% PARAMETER INPUTS
%   'parameterName' <class>: description (default: )
%
% OUTPUTS
%   out1 <class>: description
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

function [u,srtIdx] = wavetemplate(w, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'WAVETEMPLATE';

% add required, optional, and parameter-value pair arguments
addRequired(P, 'w', @iscell)
addOptional(P, 'unitNo', [], @isnumeric)
addParameter(P, 'sort', false, @islogical)
addParameter(P, 'plot', false, @islogical)

% clear workspace (parser object retains the data while staying small)
parse(P, w, varargin{:});
clear ans varargin

%% Construct template

[nUnits,nChannels] = size(w);
waveLen = size(w{1},2);

% average and reshape to length x channel x unit
u = cellfun(@(x) mean(x,1)',w,'uni',false);
u = cell2mat(reshape(u',1,nChannels,nUnits));

if P.Results.sort
    % resort by wave energy
    ener = squeeze(sum(u.^2,1));
    ener = mean(ener,1);
    [~,srtIdx] = sort(ener);
    u = u(:,:,srtIdx);
end

% Plot
if P.Results.plot
    plotwavetemplate(u)
end