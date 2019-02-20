%% Load data as Ephys object

% v: Ephys data (data points x channels)
% Fs: sample rate (Hz)

Eph = Ephys(Fs,v);

%% Inspect data to determine filter parameters

Eph.plot

chanNum = 1:2;
timeLim = [0 30]; % seconds

Eph.chan(chanNum).range(timeLim).plot('fig','clf')

Eph.chan(1).range([0 5]).plot('fig','new')

%% Filter

FILT_ORD = 2; % filter order (Butterworth)
FILT_CUT = [500 2000]; % Hz

Eph = Eph.filt('bandpass', FILT_ORD, FILT_CUT);

% Eph.filt('high',2,500);

%% Normalize

[Eph,sigma] = Eph.normalize;
% recover original data scale with Eph = Eph.data.*sigma;

%% Detect and align spikes

% initial detection
Spk = Eph.detect_spikes;

% check detection with Eph.plot('Spk',Spk);
% can adjust detection threshold with Eph.detect_spikes(thresh) with thresh
% denoting the multiplier on the standard deviation of the noise

% align
ALIGN_WIN = 3e-3; % seconds
REFRACTORY_DUR = 5e-4; % seconds

Spk = Spk.align(Eph.data, ALIGN_WIN, REFRACTORY_DUR);

% extract waveforms
WAVEFORM_DUR = 5e-3; % seconds

Spk = Spk.get_waveforms(Eph.data, WAVEFORM_DUR);

%% Cluster waveforms

NUM_FEATURES = 4;
MAX_CLUSTERS = 50;

Spk = Spk.cluster(NUM_FEATURES, MAX_CLUSTERS);

%% Inspect clustering

PLOTS_PER_PAGE = 25;
Spk(1).plot_hist('ppp',PLOTS_PER_PAGE)

%% Remove outliers

cutoff = 0.9;
Spk = Spk.triage(cutoff);

%% Compute mutual information between spike trains across channels

% ensure matched number of samples
for ii = 1:length(Spk)
    Spk(ii).nSamples = Eph.nSamples;
end

% spike indices
s = Spk.sparsify;
s = cat(2,s{:});
s = mat2cell(s,size(s,1),ones(1,size(s,2)));
spkIdx = cellfun(@find,s,'uni',false);

chanPerUnit = Spk.unit_channel;

maxLag = 5e-3; % seconds

[mi,normPkMI] = mutualinfo(spkIdx,3e4,Spk(1).nSamples,'maxLag',maxLag,'group',chanPerUnit);

%% Group clusters across channels

MI_THRESH = 0.04;
grp = findpaths(double(normPkMI > MI_THRESH));

%% Get unique sources

nSpks = cellfun(@length,spkIdx);
unitNo = zeros(length(grp),1);

for ii = 1:length(unitNo)
    [~,maxIdx] = max(nSpks(grp{ii}));
    unitNo(ii) = grp{ii}(maxIdx);
end

unqSpkIdx = spkIdx(unitNo);

%% Get wave template

WAVE_DUR = 20e-3; % seconds
waveLen = round(WAVE_DUR*Eph.Fs);

waveforms = spktrig(Eph.data,unqSpkIdx,waveLen);
w = wavetemplate(waveforms);

%% Manually prune noise templates

badUnit = [];
remUnit = setdiff(1:size(w,3), badUnit);

plotwavetemplate(w(:,:,remUnit))

% while, "dissatisfied" ---------------------------------------------------

badUnit = sort([badUnit, remUnit([2])]);
remUnit = setdiff(1:size(w,3), badUnit);
plotwavetemplate(w(:,:,remUnit))

% end ---------------------------------------------------------------------

%% Get noisy patch of data

% manual method
T_INIT = 0;
NOISE_DUR = 1; % ideal
vn = Eph.range(T_INIT+[0 NOISE_DUR-1/3e4]).data;

% automatic method
% vn = getnoise(Eph.data,Eph.Fs);

%% Estimate noise covariance

Cn = noisecov(vn,Fs,WAVE_DUR);

FUDGE_FACTOR = 1e-6; % smaller the better

Cinv = (Cn + eye(size(Cn,1))*FUDGE_FACTOR)^(-1);

%% Get final spike times

% testing
% botm(double(Eph.range([0 5]).data), Eph.Fs, w, Cinv, 'plot',{'fit','dis','res','resen'});

% w = w(:,:,remUnit);

[sEst,~,resEner] = botm(double(Eph.data),Eph.Fs,w,Cinv,'verbose',true);


