
% Adapted from ECoG_utils_private/scripts/bairECOG_convertToBIDS_TEMPLATE.m

% This script Takes BAIR data from NYU School of Medicine, gets onsets,
% writes out separate runs for each tasks, including tsv event files, and
% required BIDS metadata (coordsystem json and electrodes and channels tsv
% files). It is meant to be run cell-by-cell because some manual inputs are
% required for trigger channel selection and identifying noisy channels.

% Check whether we have the ECoG_utils repository on the path
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
task_label  = {'temporalpattern', ...            
               'temporalpattern', ... 
               'temporalpattern', ... 
               'temporalpattern', ... 
               'temporalpattern', ... 
               'temporalpattern', ... 
              };              
run_label = {'01','02','03','04','05','06'};
% NOTE: task and run labels should be noted in the order they were run!

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

% Only included the temporal task run. The fourth temporal task was a broken off run.

% Define the trigger channel name (probably a 'DC' channel, see hdr.label).
triggerChannelName = 'DC4';
triggerChannel = find(strcmp(triggerChannelName,hdr.label));
figure;plot(rawdata(triggerChannel,:)); 
title([num2str(triggerChannel) ': ' hdr.label{triggerChannel}]);
        
run_start = 720046; % Manually determined from plot of triggerchannel 
run_end   = 2019390; 

% Clip the data
data = rawdata(:,run_start:run_end);
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

% Plot the raw voltage time course of each channel
if makePlot
    for cChan = 1:1:size(data,1) 
        figure;plot(t,data(cChan,:)); 
        title([num2str(cChan) ': ' hdr.label{cChan}]);
        xlabel('Time (s)'); ylabel('Raw amplitude (microV)'); set(gca,'fontsize',16); 
        waitforbuttonpress; close; 
    end 
end
