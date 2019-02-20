%% HISTFUN histogram of functions
% Constructs a 2D histogram over an array of functions y(t). Useful for 
% visualizing the distribution of many functions, when plotting each 
% individual trace would be computationally expensive.
%
% SYNTAX
%   [counts, bins] = histfun(Y, varargin)
%
% REQUIRED INPUTS
%   Y (numeric): array of functions (observations x samples)
%
% OPTIONAL INPUTS: none
%
% PARAMETER INPUTS
%   'lim', <numeric>: length-2 vector containing the lower and upper bounds
%       on the bins. By default, uses the absolute minimum and maximum of Y
%
%   'plot', <logical>: false by default. If true, visualizes the logarithm
%       of the counts as a heatmap.
%
%   'nYTick', <integer>: number of y-axis tick marks used for plotting
%
% OUTPUTS
%   counts (numeric): 2D array of function counts at each bin and sample
%   bins (numeric): bins used for constructing the histogram
%
% EXAMPLE
%
%   % construct set of sinusoids corrupted by Gaussian noise
%   N = 10^4;
%   t = 0:1e-3:1;
%   nSin = 3;
%   amplitudes = randsample(1:25,nSin,false);
%   frequencies = randsample(1:5,nSin,false);
%   params = datasample([amplitudes;frequencies]',N,1);
%   amp = num2cell(params(:,1));
%   freq = num2cell(params(:,2));
%   Y = cell2mat(cellfun(@(A,f,sig) A*sin(2*pi*f*t)+normrnd(0,sig,size(t)), amp, freq, num2cell(randi(2,[N,1])), 'uni', false));
%
%   % plot histogram
%   histfun(Y,'plot',true);
%
%   % constrain limits and only show extrema
%   histfun(Y,'plot',true,'lim',[-30 30],'nYTick',3);
%
% IMPLEMENTATION
% Other m-files required: SMOOTH1D
% Subfunctions: SMOOTH1D
% MAT-files required: none
%
% SEE ALSO:

% Authors: Najja Marshall
% Emails: njm2149@columbia.edu
% Dated: April 2018

function [counts, bins] = histfun(Y, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'HISTFUN';

% validation functions
isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'Y', @isnumeric)
addParameter(P, 'lim', [], @(x) isempty(x) || (isnumeric(x) && length(x)==2 && x(1)<x(2)))
addParameter(P, 'plot', false, @logical)
addParameter(P, 'nYTick', 10, @(x) isscalarnum(x,0,Inf) && x==round(x))

% clear workspace (parser object retains the data while staying small)
parse(P, Y, varargin{:});
clear ans varargin

%% Histogram

% input dimensions
T = size(Y,2);

% bin limits
lim = P.Results.lim;
if isempty(lim)
    lim = double([floor(min(Y(:))), ceil(max(Y(:)))]);
end      

% bin samples
bins = lim(1):lim(2);
nBins = length(bins);

counts = cell(1,T);
for t = 1:T
    counts{t} = histc(Y(:,t),bins);
end

%% Plot

if ~P.Results.plot
    counts = cell2mat(counts);
    return
end

% filter counts
filtCounts = cell2mat(cellfun(@(x) smooth1D(x,1,'gau','sd',1), counts, 'uni', false));
counts = cell2mat(counts);

% log filtered counts
epsilon = 10^ceil(log10(abs(min(filtCounts(:)))));
logCounts = log(filtCounts/100 + epsilon);

% plot heat map
imagesc(logCounts);

% color bar
colormap('bone')
colorbar

% truncate color axis to emphasize noise
cax = caxis;
caxis([cax(1)/4, cax(2)])

% set y-axis values to increasing
set(gca,'ydir','normal')

% y-tick
firstNonNeg = find(bins>=0, 1);
yt = round(linspace(1,nBins,P.Results.nYTick));
yt = yt-(yt(find(yt>=firstNonNeg,1))-firstNonNeg);
yt = yt(yt>0 & yt<=nBins);
set(gca,'yTick',yt);

% y-tick labels
ytl = cellfun(@num2str,num2cell(round(bins(yt))),'uni',false);
set(gca,'yTickLabel',ytl);