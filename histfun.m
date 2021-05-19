%% HISTFUN histogram of functions
% Constructs a 2D histogram over an array of functions y(t). Useful for 
% visualizing the distribution of many functions, when plotting each 
% individual trace would be computationally expensive.
%
% SYNTAX
%   [counts, bins, logCounts] = histfun(Y, varargin)
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
%
%   bins (numeric): bins used for constructing the histogram
%
%   logCounts (numeric): filtered log counts used for plotting
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
%   % plot a subset of lines
%   figure
%   plot(Y(datasample(1:size(Y,1),10),:)','k')
%
%   % plot histogram of all lines
%   figure
%   histfun(Y,'plot',true);
%
%   % constrain limits and tick marks
%   figure
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

function [counts, bins, logCounts] = histfun(Y, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'HISTFUN';

% validation functions
isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'Y', @(x) isnumeric(x) || iscell(x))
addParameter(P, 'lim', [], @(x) isempty(x) || (isnumeric(x) && length(x)==2 && x(1)<x(2)))
addParameter(P, 'plot', false, @logical)
addParameter(P, 'nYTick', 5, @(x) isscalarnum(x,0,Inf) && x==round(x))
addParameter(P, 'nBins', [], @(x) isempty(x) || isnumeric(x))
addParameter(P, 'axes', [], @(x) isempty(x) || isa(x,'matlab.ui.control.UIAxes') || isa(x,'matlab.graphics.axis.Axes'))
addParameter(P, 'yTick', [], @(x) isempty(x) || isnumeric(x))
addParameter(P, 'yTickLabel', [], @(x) isempty(x) || iscell(x))
addParameter(P, 'count', [], @(x) isempty(x) || isnumeric(x))

% clear workspace (parser object retains the data while staying small)
parse(P, Y, varargin{:});
clear ans varargin

%% Histogram

if iscell(Y)
    
    Y = Y(:);
    nGrp = length(Y);
    
    % normalize
    absMax = cellfun(@(x) max(max(abs(x),[],2),[],1), Y);
    Y = cellfun(@(x,m) x/m, Y, num2cell(absMax), 'uni',false);
    
    % apply vertical offset to each group
    baseline = flipud([0;cumsum(2*ones(nGrp-1,1))]);
    Y = cellfun(@(x,b) x+b, Y, num2cell(baseline), 'uni',false);
    
    Y = cell2mat(Y);
    
    inputCell = true;
else
    inputCell = false;
end

% input dimensions
T = size(Y,2);

% bin limits
lim = P.Results.lim;
if isempty(lim)
    lim = 1*double([floor(min(Y(:))), ceil(max(Y(:)))]);
end      

% bin samples
if isempty(P.Results.nBins)
    bins = lim(1):lim(2);
    nBins = length(bins);
else
    bins = linspace(lim(1),lim(2),P.Results.nBins);
    nBins = P.Results.nBins;
end

counts = cell(1,T);
for t = 1:T
    counts{t} = histc(Y(:,t),bins);
end

%% Plot

% filter counts
filtCounts = cell2mat(cellfun(@(x) smooth1D(x,1,'gau','gauSd',1), counts, 'uni', false));
counts = cell2mat(counts);

% log filtered counts
epsilon = 10^ceil(log10(abs(min(filtCounts(:)))));
logCounts = log(filtCounts/100 + epsilon);
% logCounts = logCounts - min(logCounts(:));
% logCounts = logCounts/max(logCounts(:));


if ~P.Results.plot
    return
end

% plot heat map
if isempty(P.Results.axes)
    ax = gca;
else
    ax = P.Results.axes;
end
cla(ax,'reset')
imagesc(ax,logCounts);
set(ax,'xLim',[1 T])

% color bar
colormap(ax,'bone')
colorbar('peer',ax)

% truncate color axis to emphasize noise
cax = caxis(ax);
caxis(ax,[cax(1)/4, cax(2)])

% set y-axis values to increasing
set(ax,'ydir','normal')

% y-tick
if ~inputCell
    if isempty(P.Results.yTick)
        firstNonNeg = find(bins>=0, 1);
        yt = round(linspace(1,nBins,P.Results.nYTick));
        yt = yt-(yt(find(yt>=firstNonNeg,1))-firstNonNeg);
        yt = yt(yt>0 & yt<=nBins);
    else
        [~,yt] = min(abs(bins(:)-P.Results.yTick(:)'),[],1);
    end
    set(ax,'yLim',[1 nBins])
else
    yTickVal = 1/2;
    yTickPos = [baseline-yTickVal, baseline, baseline+fliplr(yTickVal)];
    yTick = round((yTickPos - baseline).*absMax(:));
    
    yTickPos = reshape(fliplr(yTickPos)',(1+2*length(yTickVal))*nGrp,1);
    yTick = reshape(fliplr(yTick)',(1+2*length(yTickVal))*nGrp,1);
    
    yTickPos = flipud(yTickPos);
    
    [~,yt] = min(abs(bins(:)-yTickPos(:)'),[],1);
    
    [~,yl] = min(abs(bins(:)-[-1.5 baseline(1)+1.5]),[],1);
    set(ax,'yLim',yl)
end
set(ax,'yTick',yt)

% display count
if ~isempty(P.Results.count)
    yl = get(ax,'yLim');
    text(ax,0.05*T,0.95*yl(2),sprintf('n = %i',P.Results.count),'color','w','fontsize',15)
end

% y-tick labels
if ~inputCell
    if isempty(P.Results.yTickLabel)
        ytl = cellfun(@num2str,num2cell(round(bins(yt))),'uni',false);
    else
        ytl = P.Results.yTickLabel;
    end
else
    ytl = cellfun(@num2str,flipud(num2cell(yTick)),'uni',false);
end
set(ax,'yTickLabel',ytl);

% formatting
box(ax,'off')