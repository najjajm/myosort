%% FUNCTIONTEMPLATE full function name
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

function leave1wavefit(w, Fs, xNoise, Cinv, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'LEAVE1WAVEFIT';

% validation functions
% isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'w', @isnumeric)
addRequired(P, 'Fs', @isnumeric)
addRequired(P, 'xNoise', @isnumeric)
addRequired(P, 'Cinv', @isnumeric)
% addOptional(P, 'optIn', default, validationFunction)
% addParameter(P, 'parameterName', default, validationFunction)

% clear workspace (parser object retains the data while staying small)
parse(P, w, Fs, xNoise, Cinv, varargin{:});
clear ans varargin

%%

[waveLen,nChan,nUnit] = size(w);
frameLen = 3*waveLen;
assert(size(xNoise,1) >= frameLen,'Insufficient noise data')

frameMid = round(size(xNoise,1)/2);
xn = spktrig(xNoise', {frameMid}, frameLen);
xn = cell2mat(xn')';
win = (-waveLen/2:waveLen/2-1) + round(size(xn,1)/2);

% leave-one-out template fit
spkIdx = cell(nUnit);
for ii = 1:nUnit
    
    y = xn;
    for jj = 1:nChan
        y(win,jj) = xn(win,jj) + w(:,jj,ii);
    end
    
    testUnits = setdiff(1:nUnit,ii);
    si = botm(y, Fs, w(:,:,testUnits), Cinv, 'plot', {'fit','rec','res'});
    spkIdx(ii,testUnits) = si;
    
    openFigs = findobj('Type', 'figure');
    figNo = arrayfun(@(x) x.Number,openFigs);
    figure(figNo(3))
    
    subplot(nChan,1,1)
    title(sprintf('Fit to unit %i',ii),'fontsize',15)
    subplot(nChan,1,nChan)
    if any(~cellfun(@isempty,si))
        legH = get(gca,'legend');
        parts = cellfun(@(x) strsplit(x,' '),legH.String,'uni',false);
        unitIdx = cellfun(@(x) str2double(x{2}),parts);
        legH.String = cellfun(@(x) sprintf('unit %i',testUnits(x)),num2cell(unitIdx),'uni',false);
    else
        delete(get(gca,'legend'))
    end
end