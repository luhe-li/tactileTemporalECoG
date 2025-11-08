
% % Check whether we have the ECoG_utils repository on the path
if ~exist('createBIDS_ieeg_json_nyuSOM.m')
    tbUse ECoG_utils;
end

%% Define paths and BIDS specs %%

% Input paths specs
patientID   = 726; % Specify patient's raw folder name here
RawDataDir  = '/Volumes/server/Projects/BAIR/Data/Raw/ECoG/tactile/';
BIDSDataDir = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% BIDS specs: assuming defaults for a first session, full visual set:
projectName = 'tactile';
sub_label   = 'ny726'; % Specify patient's code name here;
ses_label   = 'nyuecog01';
ses_labelt1 = 'som3t01';
acq_label   = 'clinical';
task_label  = {'tacttestascending'};              
run_label = {'01'};

% Make plots?
makePlot = 1;
% NOTE: Figures will be saved into
% derivatives/preprocessed/<sub-label>/<ses-label>/figures/bidsconversion

%% DEFINE PATHS AND DATA

% Define paths
[dataReadDir, dataWriteDir, stimWriteDir, T1WriteDir, preprocDir] = bidsconvert_getpaths(patientID, RawDataDir, ...
    BIDSDataDir, projectName, sub_label, ses_label, ses_labelt1);

% Read ECoG data
[rawdata, hdr] = bidsconvert_readecogdata(dataReadDir, ses_label);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% START OF MANUAL SECTION %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the cells below by adding several manual inputs

%% DATA TRIM (PATIENT- and SESSION-SPECIFIC)

% Only included the localizer run in the beginning.

% Define the trigger channel name (probably a 'DC' channel, see hdr.label).
triggerChannelName = 'DC4';
triggerChannel = find(strcmp(triggerChannelName,hdr.label));
figure;plot(rawdata(triggerChannel,:)); 
title([num2str(triggerChannel) ': ' hdr.label{triggerChannel}]);
        
run_start = 652639; % Manually determined from plot of trigger channel 
run_end   = 712846; 

% Clip the data
data = rawdata(:, run_start:run_end);
hdr.nSamples = size(data,2);

% Check if we have all the triggers we want
figure;plot(data(triggerChannel,:)); 
title([num2str(triggerChannel) ': ' hdr.label{triggerChannel}]);


%% BAD CHANNEL IDENTIFICATION

% Manually click through each channel to identify to trigger channel, as
% well as bad channels (specify in next cell)

% Define time axis (in seconds). First time point = 0 (this is assumed by
% the function we used to detect triggers below, and also in fieldtrip).
t = ((0:hdr.nSamples-1)/hdr.Fs); 

% % Plot the raw voltage time course of each channel
% if makePlot
%     for cChan = 1:1:size(data,1) 
%         figure;plot(t,data(cChan,:)); 
%         title([num2str(cChan) ': ' hdr.label{cChan}]);
%         xlabel('Time (s)'); ylabel('Raw amplitude (microV)'); set(gca,'fontsize',16); 
%         waitforbuttonpress; close; 
%     end 
% end

% Automatic identification
% Identify channels with values exceeding the threshold
threshold = 200;
exclude_inx = [];
for cChan = 1:size(data, 1)
    if mean(data(cChan, :)) > threshold || mean(data(cChan, :)) <- threshold || max(data(cChan,:))>1500
        exclude_inx = [exclude_inx, cChan];
    end
end

% Manually add bad channels
manual_exclude_inx = [];
exclude_inx = sort(unique([exclude_inx, manual_exclude_inx]));
% 
% % Check bad channels
for cChan = 1:1:numel(exclude_inx)
    exc_ch = exclude_inx(cChan);
    figure;plot(t,data(exc_ch,:));
    title([num2str(exc_ch) ': ' hdr.label{exc_ch}]);
    xlabel('Time (s)'); ylabel('Raw amplitude (microV)'); set(gca,'fontsize',16);
end


%% WRITE DOWN THE FOLLOWING

% Trigger channel name (probably a 'DC' channel, see hdr.label)
triggerChannelName = 'DC4'; % tactile amplifier

% Specify reasons for marked as bad, e.g. spikes, elipeptic,
% outlierspectrum, lowfreqdrift
BADCHANNELS_MANUALTABLE = [num2cell(exclude_inx)' repmat({'spikes'}, [length(exclude_inx) 1])];
badChannels = cell2mat(BADCHANNELS_MANUALTABLE(:,1));
badChannelsDescriptions = BADCHANNELS_MANUALTABLE(:,2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% END OF MANUAL SECTION %%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note that the stimuli names are recoded in the correct order, where the
% original stimulus file name is unchanged within each tsv and mat file.

% From here on, everything should run automatically:

% AUTOMATED EXTRACTION %%

% Get trigger time points from data file
peakOpts.minPeakHeight = 0.85;
peakOpts.minPeakDistance = 0.5;

[trigger_onsets] = bidsconvert_findtriggers(data, hdr, triggerChannel, makePlot);
if makePlot
    saveas(gcf, fullfile(preprocDir, 'figures', 'bidsconversion', sprintf('%s-%s-triggers_found',sub_label, ses_label)), 'epsc');
end

% Generate electrode files
[electrode_table, channel_table] = bidsconvert_getelectrodefiles(dataReadDir, hdr, triggerChannel, badChannels, badChannelsDescriptions);

% Read in stimulus files
stimDir = fullfile(dataReadDir,'localizer');
[stimData, triggersAreMatched, runTimes] = bidsconvert_matchstimulusfiles(stimDir, patientID, ses_label, task_label, run_label, trigger_onsets, 1);
if makePlot
    saveas(gcf, fullfile(preprocDir, 'figures', 'bidsconversion', sprintf('%s-%s-triggers_requested',sub_label, ses_label)), 'epsc');
end

% WRITING OF FILES %%%

% Write run files
[dataFileNames] = bidsconvert_writerunfiles(dataWriteDir, stimWriteDir, ...
    sub_label, ses_label, task_label, acq_label, run_label, ...
    data, hdr, stimData, channel_table, trigger_onsets);

% Write session files
bidsconvert_writesessionfiles(dataReadDir, dataWriteDir, T1WriteDir, ...
    sub_label, ses_label, acq_label, ses_labelt1, electrode_table, dataFileNames, runTimes);
