% This srcipt run ICA on the raw data and store results in server for later comparison

clear; close all;

restoredefaultpath
git_dir = fileparts(pwd);
ft_path = fullfile(git_dir,'fieldtrip');
addpath(ft_path)

% Call the FieldTrip defaults
ft_defaults;

%% load data

projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile';
subject           = 'ny726';
tasks             = 'temporalpattern';
sessions          = 'nyuecog01';
outputFolder      = 'ICA_reref';

% Specify data
[session, tasks, runnums] = bidsSpecifyData(projectDir, subject, sessions, tasks, []);

% Get the raw data from BIDS folder
% Data: channel x time series across 5 runs
[data, channels, events, srate] = bidsEcogGetPreprocData(projectDir, subject, sessions, tasks);

%% ICA

% save the-referenced data to the corresponding folder
outputDir = fullfile(projectDir, 'derivatives', outputFolder);
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

%% covert to fieltrip format

ieeg = [];
ieeg.label = channels.name;   % cell array of channel labels, e.g. {'V1', 'V2', ...}
ieeg.fsample = srate;         % sampling frequency
ieeg.trial = {data};          % one trial containing all channels
ieeg.time  = {(0:size(data,2)-1) ./ srate};  % time vector in seconds
ieeg.sampleinfo = [1 size(data,2)];

% write a header too
ieeg.hdr = [];
ieeg.hdr.nChans = size(channels,1);
ieeg.hdr.label  = channels.name;
ieeg.hdr.chaninfo = channels;   % keep the tables

% NOTE: events are coded 119-130; epochs are [-1 1.5];event_info = struct([]);
for i = 1:height(events)
    event_info(i).type     = events.trial_type(i);
    event_info(i).sample   = events.event_sample(i);   % sample index
    event_info(i).value    = events.trial_name{i};
    event_info(i).duration = round(events.duration(i) * ieeg.fsample);
end

fid = fopen([ 'sub-' subject '_electrodes.tsv']);
txt = textscan(fid, '%s%f%f%f%s', 'delimiter', '\t', 'Headerlines', 1);
fclose all;

% visualize data
cfg = [];
cfg.continuous = 'yes';
cfg.blocksize = 25;
cfg.preproc.detrend = 'yes';
cfg.viewmode = 'vertical';
cfg_art = ft_databrowser(cfg, ieeg);
save cfg_art cfg_art

% exclude R8!
% missing electrode positions for: 'W14', 'W15', 'H18', 'P10', 'R8' , 'S15'

%% select those channels as eeg

% save Fixed elec later 
ieeg.elec.elecpos = horzcat(txt{2:4});
ieeg.elec.label = txt{1};
cfg = [];
cfg.channel = ieeg.elec.label ;
ieeg = ft_selectdata(cfg, ieeg);

%% now sort the labels in elec correctly!!!

elec_new = ieeg.elec;

elec_new.chanpos(numel(ieeg.label)+1:end,:) = [];
elec_new.chantype(numel(ieeg.label)+1:end,:) = [];
elec_new.chanunit(numel(ieeg.label)+1:end,:) = [];
elec_new.elecpos(numel(ieeg.label)+1:end,:) = [];
elec_new.label(numel(ieeg.label)+1:end,:) = [];

for ll = 1 : numel(ieeg.label)
    ind = match_str(ieeg.elec.label, ieeg.label{ll});
    [ll, ind]
    elec_new.chanpos(ll,:) = ieeg.elec.chanpos(ind,:);
    elec_new.chantype{ll} = ieeg.elec.chantype{ind};
    elec_new.chanunit{ll} = ieeg.elec.chanunit{ind};
    elec_new.elecpos(ll,:) = ieeg.elec.elecpos(ind,:);
    elec_new.label{ll} = ieeg.elec.label{ind};
end
ieeg.elec = elec_new;

assert(...
all(cellfun(@strcmp, ieeg.elec.label, ieeg.label)), ...
'Mismatch found between ieeg.elec.label and ieeg.label');

save(fullfile(outputDir, 'ieeg.mat'), 'ieeg');

%% annotate strong artifacts that may mess up the ica computation

cfg = [];
cfg.continuous = 'yes';
cfg.blocksize = 25;
cfg.preproc.detrend = 'yes';
cfg.viewmode = 'vertical';
cfg_art = ft_databrowser(cfg, ieeg);

data_vec = ieeg.trial{1};
si_vector = [ieeg.sampleinfo(1,1):ieeg.sampleinfo(1,2)];
art_vec = zeros(size(si_vector)); % this will be a long vector holding 1 for artifactual samples and 0 for clean samples.
for aa = 1: size(cfg_art.artfctdef.visual.artifact,1)
    
    art_vec(cfg_art.artfctdef.visual.artifact(aa,1)-0.2*ieeg.fsample <= ...
        si_vector & cfg_art.artfctdef.visual.artifact(aa,2)+0.2*ieeg.fsample >= ...
        si_vector) = 1;
end

% overwrite with nan for ICA computation
data_vec(:,find(art_vec)) = nan;

%%

tst = ieeg;
tst.trial{1} = data_vec;

cfg = [];
cfg.blocksize = 25;
cfg.viewmode = 'vertical';
cfg.detrend = 'yes';
cfg.continuous = 'yes';
ft_databrowser(cfg, tst);

%% do the ica with/without!? highpass
% cfg = [];
% cfg.hpfilter    = 'yes';
% cfg.hpfreq      = 1;
% tst             = ft_preprocessing(cfg, tst);
data_comp       = ft_componentanalysis([], tst); %using eeglab runica

% pull out the unmixing matrix: 
unmixing        = data_comp.unmixing;
mixing          = inv(unmixing);
used_labels     = tst.label;

% also.. save the unmixing matrix and the used_labels
cfg = [];
cfg.unmixing = unmixing;
cfg.topolabel = ieeg.label;

% apply umixing matrix to actual data
data_comp2 = ft_componentanalysis(cfg, ieeg);

save(fullfile(outputDir, 'unmixing.mat'), 'unmixing');

% data = mixing * components 
% unmixing * data = components 

% inspect component time-courses! 
cfg          = [];
cfg.viewmode = 'vertical';
cfg.blocksize = 25;
cfg.continuous = 'yes';
ft_databrowser(cfg,data_comp2);

%% estimated the "broadness" 

chi_val = sum(bsxfun(@rdivide, ((bsxfun(@minus,abs(mixing),  mean(abs(mixing),1))).^2), ...
    mean(abs(mixing),1)),1);

%sort the chi-square values and keep track of the sorting-indices
[chisort, s_indices] = sort(chi_val, 'descend');

crit_20 = numel(find(chisort>chi2inv(0.8, numel(used_labels)-1)));

%% plot again the mixing matrix(topography) with components sorted by broadness

figure
imagesc((mixing(:, s_indices)))
caxis([-max(abs(mixing(:))) max(abs(mixing(:)))]); %axis xy
colormap(jet(256)); 
colorbar;

set(gca, 'XTick',1:numel(used_labels))
set(gca, 'YTick',1:numel(used_labels))

set(gca, 'YTickLabel', used_labels)
tmp = [1:numel(used_labels)];
set(gca, 'XTickLabel', num2cell(tmp(s_indices))); % set the xticklabels to the correct component number
colormap(jet(256));

line([crit_20+0.5, crit_20+0.5], [0, numel(used_labels)+0.5] , 'LineStyle', '--', 'LineWidth', 2, 'Color', [1 0 1]);
set(gcf, 'Color', 'w');
saveas(gcf, fullfile(outputDir, 'Sorted_Comp_LOADINGS.png'))
% remember to inspect time course of broad line components too

%% reject broad component(s)

% broad_indices = [119 120];
broad_indices = [1 120];
cfg           = [];
cfg.component = broad_indices;
ieeg_reref    = ft_rejectcomponent(cfg,data_comp2);

%% inspect re-referenced data
cfg          = [];
cfg.viewmode = 'vertical';
cfg.blocksize = 25;
cfg.continuous = 'yes';
ft_databrowser(cfg,ieeg_reref);

save(fullfile(outputDir, 'ieeg_reref.mat'), 'ieeg_reref');