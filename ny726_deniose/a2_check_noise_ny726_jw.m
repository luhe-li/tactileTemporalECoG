% tbUse tactileTemporalECoG

bidsRootPath = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile';
dataPath = fullfile(bidsRootPath, 'derivatives', 'ECoGCAR');
subject           = 'ny726';
session           = 'nyuecog01';
task              = 'temporalpattern';
[data, channels, events, srate] = bidsEcogGetPreprocData(dataPath, subject, [], task);

specs.epoch_t     = [-0.4 1.8];  % stimulus epoch window
specs.base_t      = [-0.4 -0.1]; % blank epoch window
specs.plot_ylim   = [-2 20];
% all channels on the gri
% specs.subplotdims = [4 8];
% specs.subplotidx  = 1:32;
% specs.plot_type   = 'average';
% specs.chan_names  = {'V','W','Y','Z'}; % First half of electrodes
% specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
% bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close
% out = bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); %close

chanidx1 = find(cellfun(@(x) ~isempty(x), strfind(channels.name, 'M11')));
chanidx2 = find(cellfun(@(x) ~isempty(x), strfind(channels.name, 'M12')));
chanidx3 = find(cellfun(@(x) ~isempty(x), strfind(channels.name, 'S13')));

figure, plot(data(chanidx1, :))

[epochs, t] = ecog_makeEpochs(data, events.onset, specs.epoch_t, channels.sampling_frequency(1));

% plot the mean across all epochs
figure, 
plot(t, mean(epochs(:,:,chanidx1), 2)); hold on
plot(t, mean(epochs(:,:,chanidx2), 2)); 
plot(t, mean(epochs(:,:,chanidx3), 2)); 

xline(0); yline(0)
set(gca, 'XTick', -2:1/110:1, 'XGrid','on')
xlim([-30 250]/1000)
% note that the 110 Hz oscillations starts immediately, at t=0, indicating
% electrial noise from the instrument not brain activity (unless there is a
% delay in the trigger, meaning that t=0 is actually later than the trigger
% is later than the event start by some amount

%%
% find the two-pulse trials
trial_idx = find(events.ISI > 0);
figure, plot(t, mean(epochs(:,trial_idx,chanidx1), 2));
xline(0); yline(0)


% highpass filter the to see the stimulus-locked frequency more clearly
hp1 = highpass(epochs(:,:,chanidx1),75,srate);
hp2 = highpass(epochs(:,:,chanidx2),75,srate);
hp3 = highpass(epochs(:,:,chanidx3),75,srate);
%%
figure
for ii = 1:320
    plot(t, hp1(:,trial_idx(ii)), ...
        t, hp2(:,trial_idx(ii)),...
        t, hp3(:,trial_idx(ii)), ...
        'LineWidth', 3);
    yline(0); 
    xline(0, 'LineWidth', 5)
    xline(.2, 'LineWidth', 5)
    title(ii)
    set(gca, 'XTick', -2:1/60:1, 'XGrid','on')
    xlim([-80 280]/1000)
    waitforbuttonpress
end
