%% FINDPULSES find pulses in data
% Detects pulses via peak detection on the de-noised signal. Pulses are
% detected along the largest input dimension.
%
% SYNTAX
%   s = findpulses(Y, Fs, varargin)
%
% REQUIRED INPUTS
%   Y <double array>: timeseries data. Must be 1 or 2-dimensional array. 
%   Fs < scalar numeric>: sampling frequency in Hz
%
% OPTIONAL INPUTS: NONE
%
% PARAMETER INPUTS
%   'thresh' <scalar numeric>: multiplier on estimated standard deviation 
%       of the noise (sigma). Data below thresh*sigma are set to zero
%       (default: 4)
%   'minWid' <scalar numeric>: minimum pulse width in seconds
%       (default: 1e-3)
%
% OUTPUTS
%   idx <vector>: pulse indices. If Y contains multiple channels, the
%       indices are embeded in a cell array
%   val <vector>: pulse values
%   sigma <scalar double>: estimate of noise standard deviation
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
% Dated: October 2018

function [loc,val,sigma] = findpulses(Y, Fs, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'FINDPULSES';

% validation functions
isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'Y', @(x) isnumeric(x) && isa(x,'double'))
addRequired(P, 'Fs', @(x) isscalarnum(x,0,Inf))
addOptional(P, 'thresh', 4, @(x) (length(x)==1 && isscalarnum(x,-Inf,Inf)) || (length(x)==2 && x(1)<x(2)))
addParameter(P, 'sym', true, @islogical)
addParameter(P, 'dim', 1, @(x) isscalarnum(x,0,Inf) && x==round(x))
addParameter(P, 'sigma', [], @(x) isempty(x) || isnumeric(x))
addParameter(P, 'minWid', 1e-3, @(x) isscalarnum(x,0,1))
addParameter(P, 'filt', true, @islogical)
addParameter(P, 'normalized', false, @islogical)
addParameter(P, 'OutputFormat', 'index', @(x) ismember(x,{'index','logical'}))

% clear workspace (parser object retains the data while staying small)
parse(P, Y, Fs, varargin{:});
clear ans varargin

%%

szY = size(Y);
dimY = ndims(Y);

% permute inputs to work across first dimension
shiftOrd = [P.Results.dim,setdiff(1:dimY,P.Results.dim)];
yPerm = permute(Y,shiftOrd);
szYP = size(yPerm);

% estimate standard deviation of noise
if ~P.Results.normalized
    if isempty(P.Results.sigma)
        sigma = median(abs(yPerm),1)/0.6745;
    else
        sigma = P.Results.sigma;
        szSig = size(sigma);
        assert(szSig(1)==1 && isequal(szSig(2:end),szY(2:end)));
        sigma = permute(sigma,shiftOrd);
    end
else
    sigma = ones([1,szYP(2:end)]);
end

% threshold signal
if P.Results.sym
    yRect = abs(yPerm);
    z = yRect .* (yRect > P.Results.thresh*sigma);
else
    if length(P.Results.thresh) == 1
        z = sign(P.Results.thresh)*yPerm.*(sign(P.Results.thresh)*yPerm > abs(P.Results.thresh)*sigma);
    else
        z = yPerm .* (yPerm < P.Results.thresh(1) & yPerm > P.Results.thresh(2));
    end
end

% find peaks in normalized rectified signal
z = z./sigma;
indices = cell(1,dimY-1);

[idx,val] = deal(cell([1,szYP(2:end)]));
for ii = 1:prod(szYP(2:end))
    
    [indices{:}] = ind2sub(szYP(2:end),ii);
    
    % initial pulse estimate
    [~,loc] = findpeaks(z(:,indices{:}));
    if isempty(loc)
        continue
    end
    isPeak = false(size(yPerm,1),1);
    isPeak(loc) = true;
    
    % convert rectified signal to pulse train
    s = z(:,indices{:}).*isPeak;
    
    % filter pulse train by half-minimum pulse width
    if P.Results.filt
        u = smooth1D(s,Fs,'gau',true,'sd',P.Results.minWid/2);
    else
        u = s;
    end
    
    % find pulse values and indices as peaks in train
    [~,pkLoc] = findpeaks(u,'MinPeakHeight',1);
    idx{1,indices{:}} = pkLoc;
    val{1,indices{:}} = yPerm(pkLoc,indices{:});
end

if strcmp(P.Results.OutputFormat, 'logical')
    loc = false(szYP);
    for ii = 1:prod(szYP(2:end))
        [indices{:}] = ind2sub(szYP(2:end),ii);
        loc(idx{1,indices{:}},indices{:}) = true;
    end
else
    loc = idx;
    if length(loc) == 1
        loc = loc{1};
    end
end

if length(val) == 1
    val = val{1};
end

% permute outputs to match input dimensions
shiftOrdRev = zeros(1,dimY);
for ii = 1:dimY
    shiftOrdRev(ii) = find(szYP==szY(ii));
end
sigma = permute(sigma,shiftOrdRev);
val = permute(val,shiftOrdRev);
loc = permute(loc,shiftOrdRev);