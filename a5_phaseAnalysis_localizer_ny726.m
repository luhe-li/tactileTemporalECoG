
clear;
% tbUse tactileTemporalECoG

projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile';
subject = 'ny726';

sessions = 'nyuecog01';
outputFolder = 'demean_CAR';
tasks = 'tacttestascending';
runnums = [];


% addpath(genpath('/Users/luhe/Documents/GitHub/fieldtrip/fileio/'))
% addpath(genpath('/Users/luhe/Documents/GitHub/fieldtrip/utilities/'))


% % Optinal: Check raw data spectrum before re-reference
% % [data, channels, events, ieeg_json, hdr] = bidsEcogReadFiles(projectDir, subject, sessions, tasks, '01');
% % figure;
% % spectopo(data, size(data,2),hdr.Fs);


% %% apply demean CAR (do it one)

% bidsEcogRereference(projectDir, subject, sessions, tasks, runnums, outputFolder)

%% Load data

inputFolder = 'demean_CAR';
description = 'reref';

% Load data
dataPath = fullfile(projectDir, 'derivatives', inputFolder);
[data, channels, events, srate] = bidsEcogGetPreprocData(dataPath, subject, sessions, tasks, runnums, description);

% Check power sprectrum
figure;
spectopo(data, size(data,2), srate);

%% Phase analysis

% Parameters
n_cycle = 5;
n_finger = 5;
trial_duration = events.duration(1) + events.ISI(1); 
T_cycle = trial_duration * n_finger;
f0 = 1 / T_cycle;

% Crop data
run_onset  = events.event_sample(1);
run_offset = run_onset + round(n_cycle*T_cycle*srate) - 1;
run_data   = data(:, run_onset:run_offset);
n          = size(run_data, 2);

% Detrend
detrend_data = detrend(run_data', 'linear')';

% FFT
X        = fft(run_data, [], 2);
Xpos     = X(:,1:floor(n/2));
k0       = n_cycle + 1;
amp_f0   = abs(Xpos(:,k0));

% Normalization
total_mag = sqrt(sum(abs(Xpos).^2, 2));
coh       = amp_f0 ./ total_mag;

[coh_sorted, ch_idx] = sort(coh,'descend');
topN = min(10, numel(ch_idx));
for i = 1:topN
    ci = ch_idx(i);
    fprintf('Ch-%s, coh %.3f\n', channels.name{ci}, coh(ci));
end

%%
% Inputs: detrend_data (nChan x n), srate, f0
n = size(detrend_data,2);
t = (0:n-1)/srate;
X = [sin(2*pi*f0*t') cos(2*pi*f0*t') ...
     sin(2*pi*2*f0*t') cos(2*pi*2*f0*t') ones(n,1)]; % include 2f0 if helpful
coefs = (X \ detrend_data')';                 % nChan x 5
yhat  = (coefs * X') ;                        % nChan x n
res   = detrend_data - yhat;                  
SSres = sum(res.^2,2);
SStot = sum((detrend_data - mean(detrend_data,2)).^2,2);
R2    = 1 - SSres./SStot;
[~, idx] = sort(R2, 'descend');
top = idx(1:10);
for ii = 1:10
fprintf('Top by R^2:\n'); fprintf('%s: R2=%.3f\n', channels.name{top(ii)}, R2(top(ii)));

end