% filter & normalize
[EMG,sigma] = EMG.filt('bandpass',2,[500 2000]).normalize;

% identify noise trace
[xNoise,noiseLim] = getnoise(EMG.data, EMG.Fs, 'normalized',true);

% estimate noise covariance matrix
Cn = noisecov(EMG.range(noiseLim).data, EMG.Fs, 20e-3); 
alpha = 0; % set small (e.g. 1e-8) if Cinv is poorly conditioned
Cinv = (Cn + alpha*eye(size(Cn,1)))^(-1);

% detect, align, cluster
Spk = EMG.detect_spikes...
    .align(EMG.data, 3e-3, 5e-4)...
    .get_waveforms(EMG.data, 5e-3)...
    .cluster(3, 25);

% get templates
waveforms = spktrig(EMG.data, cat(1,Spk.indices{:}), round(EMG.Fs*20e-3));

% manually select final templates using the template manager GUI
% export results as "templates"
templatemanager(EMG.data, EMG.Fs, Cinv, waveforms)

% get final spike times
load('templates')
spikes = botm(EMG.data, EMG.Fs, w, Cinv, 'verbose',true);