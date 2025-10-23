
clear; close all;
% tbUse tactileTemporalECoG;
addpath(fullfile(pwd, 'utils'));
gitDir = fileparts(pwd);
ft_path = fullfile(gitDir, 'fieldtrip');

addpath(ft_path);
ft_defaults;

projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/tactile';
subject           = 'ny726';
task              = 'temporalpattern';
session           = 'nyuecog01';
inputFolder       = 'ECoGBroadband_exclude110Hz'; 
description       = 'broadband';

figureDir = fullfile(projectDir, 'derivatives', 'ECoGFigures', 'sub-ny726','broadband_of_S1_electrodes');
if ~exist(figureDir, 'dir'), mkdir(figureDir); end

%% 1. Use Glasser atlas to select electrodes within somatosensory cortex

% Load electrodes (MNI xyz)
fid = fopen('sub-ny726_electrodes.tsv');
txt = textscan(fid, '%s%f%f%f%s', 'delimiter', '\t', 'Headerlines', 1);
fclose(fid);
elecpos = horzcat(txt{2:4});   % N × 3 numeric array [x y z] in mm

% Load Glasser atlas in MNI space
atlas = ft_read_atlas('HCP-MMP1_on_MNI152_ICBM2009a_nlin.nii.gz');
atlas = ft_convert_units(atlas, 'mm');   % make sure units are mm

% Convert electrode MNI coordinates to voxel indices
ijk = ft_warp_apply(inv(atlas.transform), elecpos); 
ijk = round(ijk);

% Clamp to valid volume range
sz = size(atlas.tissue);
ijk(:,1) = min(max(ijk(:,1),1), sz(1));
ijk(:,2) = min(max(ijk(:,2),1), sz(2));
ijk(:,3) = min(max(ijk(:,3),1), sz(3));

% Create binary mask for S1 regions
target_labels = [9 51 52]; % Glasser labels for S1 areas
voxelmask = zeros(size(atlas.tissue));
for label = target_labels
    voxelmask = voxelmask | (atlas.tissue == label);
end

% Find channels within voxelmask by sampling the mask at electrode positions
linInd = sub2ind(sz, ijk(:,1), ijk(:,2), ijk(:,3));
is_S1 = voxelmask(linInd);
elec_S1 = elecpos(is_S1, :);

% Get channel names for S1 electrodes
channel_names = txt{1};
S1_channel_names = channel_names(is_S1);

% Plot broadband timecourses for S1 electrodes
savePlot = 1;
clear specs;    
specs.epoch_t     = [-0.4 1.8]; % stimulus epoch window
specs.base_t      = [-0.4 -0.1]; % blank epoch window
specs.plot_ylim   = [-0.5, 3];

specs.plot_type   = 'averageSE';

specs.chan_names  = S1_channel_names';
specs.stim_names  = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); close
specs.stim_names  = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); close


%% 2. Use the anatomy label to select postcentral electrode contacts

filePath = 'autoNY726_2_T1_depth_split_STG_MTG_AnatomicalRegions.txt';

% Read all non-commented lines
rawLines = readlines(filePath);
rawLines = rawLines(~startsWith(rawLines, "%") & strlength(rawLines) > 0);

% Preallocate output cell array
nLines = numel(rawLines);
T = table('Size', [nLines, 5], ...
    'VariableTypes', {'string', 'double', 'double', 'double', 'cell'}, ...
    'VariableNames', {'Electrode', 'X', 'Y', 'Z', 'Regions'});

% Parse each line
for i = 1:nLines
    line = strtrim(rawLines(i));
    tokens = regexp(line, '\s+', 'split');
    
    % first 4 columns: Electrode, X, Y, Z
    T.Electrode(i) = string(tokens{1});
    T.X(i) = str2double(tokens{2});
    T.Y(i) = str2double(tokens{3});
    T.Z(i) = str2double(tokens{4});
    
    % The rest: alternating percentage + label pairs
    remaining = tokens(5:end);
    nPairs = floor(numel(remaining)/2);
    
    regionList = cell(nPairs, 2);
    for j = 1:nPairs
        pctStr = remaining{2*j-1};
        labelStr = remaining{2*j};
        
        % Remove '%' and convert to number
        pctNum = str2double(erase(pctStr, '%'));
        
        regionList{j,1} = pctNum;
        regionList{j,2} = labelStr;
    end
    T.Regions{i} = regionList; % store as {percentage, label}
end

% Initialize arrays to store channel names and anatomical labels
anat_channels = {};
anat_labels = {};

% Loop through each row of the table
for i = 1:height(T)
    regions = T.Regions{i};
    for j = 1:size(regions,1)
        if contains(regions{j,2}, 'postcentral', 'IgnoreCase', true)
            % Add this electrode name to our list
            electrode = T.Electrode(i);
            if regexp(electrode, '[A-Z]0\d')
                electrode = regexprep(electrode, '0', '');
            end
            anat_channels{end+1} = electrode;
            anat_labels{end+1} = regions{j,2};
            break;
        end
    end
end

% Plot broadband timecourses for postcentral electrodes
specs.chan_names  = cellfun(@char, anat_channels, 'UniformOutput', false); % Convert to char array format
savePlot = 1;

specs.stim_names = {'ONE-PULSE-1', 'ONE-PULSE-2', 'ONE-PULSE-3', 'ONE-PULSE-4', 'ONE-PULSE-5', 'ONE-PULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); close
specs.stim_names = {'TWO-PULSE-1', 'TWO-PULSE-2', 'TWO-PULSE-3', 'TWO-PULSE-4', 'TWO-PULSE-5', 'TWO-PULSE-6'};
bidsEcogPlotTrials(projectDir, subject, session, task, [], inputFolder, description, specs, savePlot); close
