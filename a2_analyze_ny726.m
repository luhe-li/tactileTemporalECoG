
clear;
tbUse tactileTemporalECoG

projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile';
subject = 'ny726';

sessions = 'nyuecog01';
outputFolder = 'demean_CAR';
tasks = 'tacttestascending';
runnums = [];

% addpath(genpath('/Users/luhe/Documents/GitHub/fieldtrip/fileio/'))
% addpath(genpath('/Users/luhe/Documents/GitHub/fieldtrip/utilities/'))

%% apply demean CAR (do it one)

bidsEcogRereference(projectDir, subject, sessions, tasks, runnums, outputFolder)

%% 
nCycle = 5; % Number of repeats across this run
