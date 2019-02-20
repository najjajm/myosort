%% FINDPULSES find pulses
% Detects pulse events in time series data along any singular dimension. 
% First, an estimate of the standard deviation of the noise [Donoho & 
% Johnstone, 1994] is used to determine a threshold for denoising the data.
% Then, a sequence of peak finding and smoothing operations are used to 
% identify the occurence of pulses that can have mono- or multiphasic 
% signatures.
%
% SYNTAX
%   [loc,val,sigma] = findpulses(X, Fs, [thresh], varargin)
%
% REQUIRED INPUTS
%   X (double): time series array
%   Fs (scalar): sample frequency in Hz
%
% OPTIONAL INPUTS:
%   thresh (numeric): threshold (in multiples of the estimated standard
%       deviation of the noise) that is used to denoise the data. Can
%       either specify a single or double threshold.
%
% PARAMETER INPUTS
%   'sym', <logical>: if true (default), the polarity of the data and the
%       threshold are ignored for denoising. Otherwise, if the threshold is
%       negative, data above it are ignored for pulse detection; if the
%       threshold is positive, data below it are ignored.
%
%   'dim', <integer>: dimension along which pulses are detected
%
%   'sigma', <numeric>: standard deviation of the noise. If empty
%       (default), this is estimated automatically. Otherwise, this can be
%       provided if known in advance, but must match the size of the data
%       array when collapsed along the detection dimension ("dim").
%
%   'minWid', <scalar>: minimum width of the pulse signature in seconds.
%       This is used to smooth the denoised signal in order to account for 
%       multiphasic features in the pulse signature. If set to 0, any peak
%       in the data will be detected as a distinctive pulse (default: 1e-3)
%
%   'normalized', <logical>: indicates whether the data has been normalized
%       by the standard deviation of the noise. If false (default),
%       estimates or uses the provided standard deviation of the noise to 
%       denoise the signal. Otherwise, assumes the standard deviation is 
%       unity for all channels.
%
%   'outputFormat', <string>: if 'index' (default), retuns the pulse
%       locations as indices. If 'logical', returns a boolean array of the
%       same dimensionality as the input with true values at each pulse
%       location.
%   
% OUTPUTS
%   loc (numeric or cell): indices of pulse locations. If X is a vector,
%       loc will be returned as a numeric vector. If X contains multiple 
%       channels, the pulse locations are embedded in a cell array.
%
%   val (numeric or cell): values of the data at each pulse index. Data
%       type depends on the dimensionality of X, as described for loc.
%
%   sigma (numeric): estimate of noise standard deviation
%
% EXAMPLES
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

function [loc,val,sigma] = findpulses(X, Fs, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'FINDPULSES';

% validation functions
isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'X', @(x) isnumeric(x) && isa(x,'double'))
addRequired(P, 'Fs', @(x) isscalarnum(x,0,Inf))
addOptional(P, 'thresh', 4, @(x) (length(x)==1 && isscalarnum(x,-Inf,Inf)) || (length(x)==2 && x(1)<x(2)))
addParameter(P, 'sym', true, @islogical)
addParameter(P, 'dim', 1, @(x) isscalarnum(x,0,ndims(X)+1) && x==round(x))
addParameter(P, 'sigma', [], @(x) isempty(x) || isnumeric(x))
addParameter(P, 'minWid', 1e-3, @(x) isscalarnum(x,-eps,Inf))
addParameter(P, 'normalized', false, @islogical)
addParameter(P, 'outputFormat', 'index', @(x) ismember(x,{'index','logical'}))

% clear workspace (parser object retains the data while staying small)
parse(P, X, Fs, varargin{:});
clear ans varargin

%% Process data

szX = size(X);
dimX = ndims(X);

% permute inputs to work across first dimension
shiftOrd = [P.Results.dim,setdiff(1:dimX,P.Results.dim)];
xPerm = permute(X,shiftOrd);
szXP = size(xPerm);

% estimate standard deviation of noise
if ~P.Results.normalized
    if isempty(P.Results.sigma)
        sigma = median(abs(xPerm),1)/0.6745;
    else
        sigma = P.Results.sigma;
        szSig = size(sigma);
        assert(szSig(1)==1 && isequal(szSig(2:end),szX(2:end)),'Mismatch between input data and noise standard deviation');
        sigma = permute(sigma,shiftOrd);
    end
else
    sigma = ones([1,szXP(2:end)]);
end

% threshold signal
if P.Results.sym
    xRect = abs(xPerm);
    z = xRect .* (xRect > abs(P.Results.thresh(1))*sigma);
else
    if length(P.Results.thresh) == 1
        z = sign(P.Results.thresh)*xPerm.*(sign(P.Results.thresh)*xPerm > abs(P.Results.thresh)*sigma);
    else
        z = xPerm .* (xPerm < P.Results.thresh(1) & xPerm > P.Results.thresh(2));
    end
end

%% Find pulses

z = z./sigma;
indices = cell(1,dimX-1);

[idx,val] = deal(cell([1,szXP(2:end)]));
for ii = 1:prod(szXP(2:end))
    
    [indices{:}] = ind2sub(szXP(2:end),ii);
    
    % initial pulse estimate
    [~,loc] = findpeaks(z(:,indices{:}));
    if isempty(loc)
        continue
    end
    isPeak = false(size(xPerm,1),1);
    isPeak(loc) = true;
    
    % convert rectified signal to pulse train
    s = z(:,indices{:}).*isPeak;
    
    % filter pulse train by half-minimum pulse width
    if P.Results.minWid > 0
        u = smooth1D(s,Fs,'gau',true,'sd',P.Results.minWid/2);
    else
        u = s;
    end
    
    % find pulse values and indices as peaks in train
    [~,pkLoc] = findpeaks(u,'MinPeakHeight',1);
    idx{1,indices{:}} = pkLoc;
    val{1,indices{:}} = xPerm(pkLoc,indices{:});
end

%% Post-processing

% format outputs
if strcmp(P.Results.outputFormat, 'logical')
    loc = false(szXP);
    for ii = 1:prod(szXP(2:end))
        [indices{:}] = ind2sub(szXP(2:end),ii);
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
shiftOrdRev = zeros(1,dimX);
for ii = 1:dimX
    shiftOrdRev(ii) = find(szXP==szX(ii));
end
sigma = permute(sigma,shiftOrdRev);
val = permute(val,shiftOrdRev);
loc = permute(loc,shiftOrdRev);