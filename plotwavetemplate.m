%% PLOTWAVETEMPLATE plot waveform template
% Plots a set of waveform templates in a grid of size channels x units.
%
% SYNTAX
%   plotwavetemplate(u, varargin)
%
% REQUIRED INPUTS
%   u (numeric): template of dimensions [wave length x channels x units]
%
% OPTIONAL INPUTS: none
%
% PARAMETER INPUTS
%   'cmap', <char>: brewermap-valid color map (default: 'Spectral')
%
%   'bright', <scalar>: constant offset (between +/- 1) applied to color 
%       brightness. Values closer to -1 darken the color map, whereas
%       values closer to +1 brighten the color map. (default: -0.2)
%
%   'xPad', <scalar>: length of padding between each column in proportion
%       of wavelength (default: 0.1)
%
%   'yPad', <scalar>: length of padding between each row in proportion
%       of the amplitude range of the largest waveform (default: 0.05)
%
%   'yTick', <numeric>: non-negative, monotonically decreasing vector that
%       specifies the y-tick locations for the templates on each channel in
%       units of fractions of the maximum template value on each channel
%       (default: [2/3 1/3])
%
%   'chanNo', <numeric>: vector of channel numbers used for text labels
%
% OUTPUTS: none
%
% EXAMPLE(S) 
%
%
% IMPLEMENTATION
% Other m-files required: BREWERMAP
% Subfunctions: BREWERMAP
% MAT-files required: none
%
% SEE ALSO: WAVETEMPLATE, BREWERMAP

% Authors: Najja Marshall
% Emails: njm2149@columbia.edu
% Dated: February 2019

function plotwavetemplate(u, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'PLOTWAVETEMPLATE';

% validation functions
isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'u', @isnumeric)
addParameter(P, 'cmap', 'Spectral', @(x) ischar(x) || (isnumeric(x) && size(x,2)==3))
addParameter(P, 'bright', -0.2, @(x) isscalarnum(x,-1,1))
addParameter(P, 'xPad', 0.1, @(x) isscalarnum(x,0,Inf))
addParameter(P, 'yPad', 0.1, @(x) isscalarnum(x,0,Inf))
addParameter(P, 'yTick', [2/3 1/3], @(x) isnumeric(x) && all(x<=1) && all(x>0) && all(diff(x)<0))
addParameter(P, 'chanNo', 1:size(u,2), @isnumeric)
addParameter(P, 'xTickLabel', [], @iscell)
addParameter(P, 'chanLabels', 'on', @(x) ischar(x) && ismember(x,{'on','off'}))
addParameter(P, 'axes', [], @(x) isempty(x) || isa(x,'matlab.ui.control.UIAxes') || isa(x,'matlab.graphics.axis.Axes'))

% clear workspace (parser object retains the data while staying small)
parse(P, u, varargin{:});
clear ans varargin

%%

[waveLen,nChan,nUnit] = size(u);

% normalize templates
if all(u==0)
    yAbsLim = ones(1,nChan);
else
    yAbsLim = max(max(abs(u),[],1),[],3);
    u = u./yAbsLim;
end

% horizontal and vertical offsets
xPos = [0,cumsum(repmat((1+P.Results.xPad)*size(u,1),1,nUnit-1))];
yPos = flipud([0;cumsum((1+P.Results.yPad)*2*ones(nChan-1,1))]);

% y-tick position and values
yTickVal = P.Results.yTick;
yTickPos = [yPos-yTickVal, yPos, yPos+fliplr(yTickVal)];
yTick = round((yTickPos - yPos).*yAbsLim(:));

yTickPos = reshape(fliplr(yTickPos)',(1+2*length(yTickVal))*nChan,1);
yTick = reshape(fliplr(yTick)',(1+2*length(yTickVal))*nChan,1);

% color map
cmap = P.Results.cmap;
if ischar(cmap)
    cmap = brewermap(nUnit,cmap);
    cmap = min([1 1 1],max([0 0 0], cmap + P.Results.bright));
end

% plot shapes
if isempty(P.Results.axes)
    figure
    ax = gca;
else
    ax = P.Results.axes;
end
cla(ax,'reset')
hold(ax,'on')
for ch = 1:nChan
    for un = 1:nUnit
        plot(ax,xPos(un)+(1:waveLen),u(:,ch,un)+yPos(ch),'color',cmap(un,:));
    end
end

% x labels
set(ax,'xtick',xPos+waveLen/2)
unitNo = num2cell(1:nUnit);
if isempty(P.Results.xTickLabel)
    set(ax,'xticklabel',cellfun(@(x) sprintf('unit %i',x),unitNo,'uni',false))
else
    set(ax,'xticklabel',P.Results.xTickLabel)
end

% y labels
set(ax,'ytick',flipud(yTickPos))
set(ax,'yTickLabel',cellfun(@num2str,flipud(num2cell(yTick)),'uni',false))

% channel numbers
if strcmp(P.Results.chanLabels,'on')
    for ch = 1:nChan
        text(ax,1.01*(waveLen+xPos(end)),yPos(ch),['channel ' num2str(P.Results.chanNo(ch))],'fontsize',14)
    end
end

set(ax,'fontsize',14)

axis(ax,'tight')
set(ax,'xlim',[0, xPos(end)+waveLen])
set(ax,'ylim',[-1, yPos(1)+1])