%% MYOSORT wrapper for myosort toolbox
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

function sort(X, Fs, varargin)

persistent RunTimer

%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'MYO.SORT';

% validation functions
% isscalarnum = @(x,lb,ub) isscalar(x) && isnumeric(x) && x>lb && x<ub;

% add required, optional, and parameter-value pair arguments
addRequired(P, 'X', @isnumeric)
addRequired(P, 'Fs', @isscalar)
% setup
addParameter(P, 'block', ':', @ischar)
addParameter(P, 'savePath', '', @ischar)
addParameter(P, 'visualize', false, @islogical)
addParameter(P, 'maxArray', 1e9, @isscalar)
% pre-processing
addParameter(P, 'filterOrder', 2, @isscalar)
addParameter(P, 'filterLowCut', 500, @isnumeric)
addParameter(P, 'filterHighCut', Inf, @isnumeric)
% detection
addParameter(P, 'detectThreshold', 6, @isscalar)
% alignment
addParameter(P, 'waveformDuration', [3, 10]*1e-3, @isnumeric)
addParameter(P, 'alignBatchSize', 100, @isscalar)
% clustering
addParameter(P, 'clusterFeatures', 3, @isnumeric)
% triage
addParameter(P, 'triageCutProportion', 0.2, @isnumeric)
addParameter(P, 'triageBins', 100, @isscalar)
addParameter(P, 'triageMinCount', 5, @isnumeric)
% merge
addParameter(P, 'mergeThreshold', 0.7, @isscalar)
% noise covariance
addParameter(P, 'noiseAlpha', 0.5, @isscalar)
% vetting
addParameter(P, 'vetCVThreshold', 1, @isscalar)
addParameter(P, 'vetMisalignThreshold', 0.1, @isnumeric)
addParameter(P, 'vetResidualThreshold', 1.1, @isscalar)
addParameter(P, 'vetBatchDuration', 10, @isscalar)
addParameter(P, 'vetRounds', 30, @isscalar)
% template matching
addParameter(P, 'refractoryDuration', 4e-3, @isscalar)
% visualization
addParameter(P, 'maxSubplotsPerPage', 5, @isscalar)

% clear workspace (parser object retains the data while staying small)
parse(P, X, Fs, varargin{:});
clear ans varargin

%% Setup

% copy input parser to sort options
opt = P.Results;
opt = rmfield(opt,'X');
clear P

% sort blocks and status updates
BlockStatus = struct(...,
    'process',          'Pre-processing data',...
    'detect',           'Detecting spikes',...
    'align',            'Aligning waveforms',...
    'cluster',          'Clustering',...
    'merge',             'Merging templates',...
    'triage',           'Triaging outliers',...
    'batch',            'Splitting data into batches',...
    'noisecov',         'Estimating noise covariance',...
    'vet',              'Vetting templates',...
    'templatematching', 'Extracting final spike times');

% timer
sortBlocks = fieldnames(BlockStatus);
RunTimer = cell2struct(cell(1+length(sortBlocks),1),[sortBlocks; 'total'],1);
RunTimer.total = tic;

% cross reference specified sort range

sortRange = strsplit(opt.block,':');
if isempty(sortRange{1})
    sortRange{1} = sortBlocks{1};
end
if isempty(sortRange{end})
    sortRange{2} = sortBlocks{end};
end
for ii = 1:2    
    assert(ismember(sortRange{2},sortBlocks), 'Unrecognized sort block: %s', sortRange{ii})
end
assert(find(ismember(sortBlocks,sortRange{1})) <= find(ismember(sortBlocks,sortRange{2})),...
    'Block %i comes before %i!', sortRange{1}, sortRange{2})

% ensure preceeding blocks have data

blockStart = find(ismember(sortBlocks,sortRange{1}));
blockStop = find(ismember(sortBlocks,sortRange{2}));

% force horizontal
if size(X,1) > size(X,2)
    warning('Transposing data array')
    X = X';
end

% ensure even wavelengths
waveLen = ceil(opt.waveformDuration*Fs);
opt.waveformDuration = (waveLen + mod(waveLen,2))/Fs;

% ensure noise covariance matrix exists
noiseCovIdx = find(ismember(sortBlocks,'noise_cov'));
if blockStart > noiseCovIdx
    if ~exist([opt.savePath 'noise_cov.mat'],'file')
        warning('Missing noise covariance matrix')
        blockStart = noiseCovIdx;
    end
end

% set save path for processed results
savePath = opt.savePath;
if isempty(savePath)
    str = dbstack;
    selfName = which(['myo.' str.name]);
    fileParts = strsplit(selfName,filesep);
    selfPath = cell2mat(join(fileParts(1:end-1),filesep));
    savePath = [selfPath, filesep, 'myosort-out', filesep];
end
opt.savePath = savePath;
if ~exist(opt.savePath,'dir')
    mkdir(opt.savePath)
end

% save sort options
save([opt.savePath 'settings'],'opt')

%% Main

for iBlock = blockStart:blockStop
    
    RunTimer.(sortBlocks{iBlock}) = tic;
    
    % open summary
    sumPref = [opt.savePath 'summary_' sortBlocks{iBlock}];
    diary([sumPref '_working.txt'])
    
    % print status
    fprintf('%s . . .\n', BlockStatus.(sortBlocks{iBlock}))
    
    % run block
    [X,opt] = myo.(sortBlocks{iBlock})(X,opt);
    
    % log runtime
    runTime = toc(RunTimer.(sortBlocks{iBlock}));
    if runTime < 60
        tUnit = 'sec';
    elseif runTime >= 60 && runTime < 3600
        runTime = runTime/60;
        tUnit = 'min';
    else
        runTime = runTime/3600;
        tUnit = 'hr';
    end
    fprintf('runtime: %.2f %s\n\n',runTime,tUnit);
    
    % close summary
    diary off
    movefile([sumPref '_working.txt'],[sumPref '.txt'])
end

%% Post-processing

% concatenate summary files
summaryFiles = dir([opt.savePath 'summary_*']);
summaryName = arrayfun(@(x) x.name,summaryFiles,'uni',false);
summaryPart = cellfun(@(x) strsplit(x,'_'),summaryName,'uni',false);
summaryPart = cellfun(@(x) x{2}, summaryPart, 'uni',false);
summaryPart = cellfun(@(x) strsplit(x,'.'),summaryPart,'uni',false);
summaryPart = cellfun(@(x) x{1}, summaryPart, 'uni',false);

writeOrder = zeros(size(summaryPart));
for ii = 1:length(summaryPart)
    writeOrder(ii) = find(ismember(sortBlocks,summaryPart{ii}));
end
[~,srtIdx] = sort(writeOrder);
summaryFiles = summaryFiles(srtIdx);

if exist([opt.savePath 'summary.txt'],'file')
    delete([opt.savePath 'summary.txt'])
end
fid = fopen([opt.savePath 'summary.txt'],'a');
for ii = 1:length(summaryFiles)
    fidPart = fopen([opt.savePath, summaryFiles(ii).name]);
    fwrite(fid,fread(fidPart,'*char'),'*char');
end

% log runtime
runTime = toc(RunTimer.total);
if runTime < 60
    tUnit = 'sec';
elseif runTime >= 60 && runTime < 3600
    runTime = runTime/60;
    tUnit = 'min';
else
    runTime = runTime/3600;
    tUnit = 'hr';
end
fwrite(fid,sprintf('Total runtime: %.2f %s',runTime,tUnit),'*char');

fclose(fid);