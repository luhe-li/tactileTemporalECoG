clear; close all;

restoredefaultpath
git_dir = fileparts(pwd);
ft_path = fullfile(git_dir,'fieldtrip');
addpath(ft_path)

% Call the FieldTrip defaults
ft_defaults;

%% load 
load('ny726_data');

%% double check, just do an fft of the raw data

[n_channel, n_sample] = size(data);
nfft = 2^nextpow2(n_sample);
freqs = (0:nfft-1) * (srate/nfft);
freq_idx = freqs >= 1 & freqs <= 200;  % select 1–200 Hz range

% Compute FFT
fft_data = fft(data, nfft, 2);
amp_data = abs(fft_data/n_sample);
amp_data = amp_data(:, 1:nfft/2+1);

figure;
plot(freqs(freq_idx), amp_data(:, freq_idx));
xlabel('Frequency (Hz)');
ylabel('Power Spectrum');

%% covert to fieltrip format

ieeg = [];
ieeg.label = channels.name;   % cell array of channel labels, e.g. {'V1', 'V2', ...}
ieeg.fsample = 512;           % sampling frequency
ieeg.trial = {data};          % one trial containing all channels
ieeg.time  = {(0:size(data,2)-1) ./ srate};  % time vector in seconds
ieeg.sampleinfo = [1 size(data,2)];

% write a header too
ieeg.hdr = [];
ieeg.hdr.nChans = size(channels,1);
ieeg.hdr.label  = channels.name;
ieeg.hdr.chaninfo = channels;   % keep the tables
% NOTE: events are coded 119-130; epochs are [-1 1.5];

% Convert event table into FieldTrip-style events
event_info = struct([]);
for i = 1:height(events)
    event_info(i).type     = events.trial_type(i);
    event_info(i).sample   = events.event_sample(i);   % sample index
    event_info(i).value    = events.trial_name{i};
    event_info(i).duration = round(events.duration(i) * ieeg.fsample);
end

data_vec = ieeg.trial{1};
tst = ieeg;
tst.trial{1} = data_vec;

%% load ica component analysis from sebastian

comp_dir = fullfile(pwd, 'components');
load(fullfile(comp_dir, 'data_comp.mat'));
load(fullfile(comp_dir, 'unmixing.mat'));

%% do the remaining

% pull out the unmixing matrix: 
mixing          = inv(unmixing);
used_labels     = tst.label;

% also.. save the unmixing matrix and the used_labels
cfg = [];
cfg.unmixing = unmixing;
cfg.topolabel = ieeg.label;

% apply umixing matrix to actual data
data_comp2 = ft_componentanalysis(cfg, ieeg);
