%% PLOTWAVETEMPLATE plot waveform template
% Plots a set of waveform templates in a grid of size channels x units
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
addParameter(P, 'cmap', 'Spectral', @ischar)
addParameter(P, 'bright', -0.2, @(x) isscalarnum(x,-1,1))
addParameter(P, 'xPad', 0.1, @(x) isscalarnum(x,0,Inf))
addParameter(P, 'yPad', 0.05, @(x) isscalarnum(x,0,Inf))

% clear workspace (parser object retains the data while staying small)
parse(P, u, varargin{:});
clear ans varargin

%%

[waveLen,nChan,nUnit] = size(u);

% color map
cmap = brewermap(nUnit,P.Results.cmap);
cmap = min([1 1 1],max([0 0 0], cmap + P.Results.bright));

% horizontal offset
xPos = [0,cumsum(repmat((1+P.Results.xPad)*size(u,1),1,nUnit-1))];

% vertical offset
yLim = [floor(min(u(:))), ceil(max(u(:)))];
yPos = round(flipud([0;cumsum(repmat((1+P.Results.yPad)*diff(yLim),nChan-1,1))]));

% y-tick position and values
yTickPos = round([yPos+yLim(2)/2, yPos, yPos+yLim(1)/2]);
yTick = yTickPos - yPos;

yTickPos = reshape(yTickPos',3*nChan,1);
yTick = reshape(yTick',3*nChan,1);

figure
hold on
for ch = 1:nChan
    for un = 1:nUnit
        plot(xPos(un)+(1:waveLen),u(:,ch,un)+yPos(ch),'color',cmap(un,:));
    end
end

% x labels
set(gca,'xtick',xPos+waveLen/2)
unitNo = num2cell(1:nUnit);
set(gca,'xticklabel',cellfun(@(x) sprintf('unit %i',x),unitNo,'uni',false))

% y labels
set(gca,'ytick',flipud(yTickPos))
set(gca,'yTickLabel',cellfun(@num2str,flipud(num2cell(yTick)),'uni',false))

% channel numbers
for ch = 1:nChan
    text(1.01*(waveLen+xPos(end)),yPos(ch),['channel ' num2str(ch)],'fontsize',14)
end

set(gca,'fontsize',14)

axis tight
set(gca,'xlim',[0, xPos(end)+waveLen])
set(gca,'ylim',[yLim(1), yPos(1)+yLim(2)])