%% SPKTRIG spike-triggered extraction
% Extracts data values aligned to a set of spike indices. Any spikes for
% which the extraction window is outside the range of the data are removed
% and the final set of spike indices are also returned.
%
% SYNTAX
%   [y, s] = spktrig(X, s, len, varargin)
%
% REQUIRED INPUTS
%   X (numeric): time series array of dimensions [channels x samples]
%   s (cell): array of spike indices
%   len (integer): length of data frame aligned about each spike
%
% OPTIONAL INPUTS: none
%
% PARAMETER INPUTS
%   'align', <char>: alignment point of the data frame relative to each
%       spike. Options: 'left', 'right', or 'center' (default).
%
% OUTPUTS
%   y (cell): array of extracted data snippets
%   s (cell): array of final spike indices
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

function [y, s] = spktrig(X, s, len, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'SPKTRIG';

% validation functions
isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'X', @isnumeric)
addRequired(P, 's', @iscell)
addRequired(P, 'len', @(x) isscalarnum(x,0,Inf) && x==round(x))
% addOptional(P, 'optIn', default, validationFunction)
addParameter(P, 'align', 'center', @(x) ismember(x,{'center','left','right'}))

% clear workspace (parser object retains the data while staying small)
parse(P, X, s, len, varargin{:});
clear ans varargin

%%

% force horizontal
if size(X,1) > size(X,2)
    X = X';
end
[nChannels, nSamples] = size(X);

% ensure length is even
len = len + mod(len,2);

% alignment
switch P.Results.align
    case 'left'
        frame = 0:len-1;
        spkLim = [1, nSamples-len+1];
        
    case 'center'
        frame = -len/2:len/2-1;
        spkLim = [1+len/2, nSamples-len/2+1];
        
    case 'right'
        frame = -len:-1;
        spkLim = [len, nSamples+1];
end

% format spike indices
s = cellfun(@(z) z(:),s,'uni',false);

% remove unbounded spike indices
s = cellfun(@(z) z(z>=spkLim(1) & z<=spkLim(2)), s, 'uni', false);
s(cellfun(@isempty,s)) = [];
nObs = length(s);

% extract data at each spike
y = cell(nObs,nChannels);
for un = 1:nObs
    for ch = 1:nChannels
        y{un,ch} = cell2mat(cellfun(@(idx) X(ch,idx+frame), num2cell(s{un}),'uni',false));
    end
end

% extract singular output from cell
if nObs == 1 && nChannels == 1
    y = y{:};
end
