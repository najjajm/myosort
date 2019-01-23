%% STA spike-triggered average
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

function [y,spkIdxNew] = sta(X, s, len, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'STA';

% validation functions
isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'X', @isnumeric)
addRequired(P, 's', @(z) iscell(z) | (islogical(z) && isequal(size(X),size(z))))
addRequired(P, 'len', @(x) isscalarnum(x,0,Inf))
% addOptional(P, 'optIn', default, validationFunction)
addParameter(P, 'align', 'center', @(x) ismember(x,{'center','left','right'}))

% clear workspace (parser object retains the data while staying small)
parse(P, X, s, len, varargin{:});
clear ans varargin

%%

% force horizontal
if size(X,1) > size(X,2)
    X = X';
    if islogical(s)
        s = s';
    end
end
[nChannels, nDataPoints] = size(X);

% alignment
switch P.Results.align
    case 'left'
        frame = 0:len-1;
        spkLim = [1, nDataPoints-len+1];
        
    case 'center'
        frame = -len/2:len/2-1;
        spkLim = [1+len/2, nDataPoints-len/2+1];
        
    case 'right'
        frame = -len:-1;
        spkLim = [len, nDataPoints+1];
end

% read spikes as indices
if iscell(s)
    spkIdx = s;
else
    spkIdx = cellfun(@find, mat2cell(s,ones(nChannels,1),nDataPoints),'uni',false);
end
spkIdx = cellfun(@(z) z(:),spkIdx,'uni',false);
nUnits = length(spkIdx);

% bound spike indices within limits
spkIdxNew = cellfun(@(z) z(z>=spkLim(1) & z<=spkLim(2)), spkIdx,'uni',false);

% take spike-triggered average
y = cell(nUnits,nChannels);
for un = 1:nUnits
    for ch = 1:nChannels
        y{un,ch} = cell2mat(cellfun(@(idx) X(ch,idx+frame), num2cell(spkIdxNew{un}),'uni',false));
    end
end

% extract singular output from cell
if nUnits == 1 && nChannels == 1
    y = y{:};
end
