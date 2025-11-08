
% % Check whether we have the ECoG_utils repository on the path
% if ~exist('createBIDS_ieeg_json_nyuSOM.m')
%     tbUse ECoG_utils;
% end

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
task_label  = {'loaclizer'};              
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

% Only included the temporal task run. Excluded the broken-off 4th run.

% Define the trigger channel name (probably a 'DC' channel, see hdr.label).
triggerChannelName = 'DC4';
triggerChannel = find(strcmp(triggerChannelName,hdr.label));
figure;plot(rawdata(triggerChannel,:)); 
title([num2str(triggerChannel) ': ' hdr.label{triggerChannel}]);
        
run_start = 736627; % Manually determined from plot of trigger channel 
t2 = 1327930;
t3 = 1414830;
run_end   = 1993870; 

% Clip the data
clip1 = rawdata(:,run_start:t2);
clip2 = rawdata(:,t3:run_end);
data = [clip1, clip2];
hdr.nSamples = size(data,2);

% Check if we have all the triggers we want
figure;plot(data(triggerChannel,:)); 
title([num2str(triggerChannel) ': ' hdr.label{triggerChannel}]);