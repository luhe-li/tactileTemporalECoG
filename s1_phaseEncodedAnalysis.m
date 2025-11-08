% clear the workspace
clear; close all; clc;

% add path to freesurfer directory
restoredefaultpath
addpath(genpath('/Applications/freesurfer'))
addpath(genpath('~/matlab/toolboxes/jsonlab')); % savejson
addpath(genpath('~/matlab/toolboxes/MRI_tools')); % bidsSpecifyEPIs

%% Setup params
projectName     = 'TemporalTactileCounting';
taskNames       = {'loc'};
dataFolder      = 'fmriprep';
dataStr         = 'space-T1w*preproc_bold.nii.gz';
allDesigns      = {'Sweep'};
whichDesign     = 1;
fingerNames     = {'Thumb', 'Index', 'Middle', 'Ring', 'Little'};
numFingers      = numel(fingerNames);
recompute       = 0;

stimdur         = 5; % seconds
nCycles         = 10;
durCycle        = stimdur * numFingers;
tr              = 1;
dropFrames      = 5;
shift           = 3;
durScan         = durCycle * nCycles + (2*dropFrames);


allModalities   = {'tactile'};
whichModality   = 1;
whichAnalysis   = 'Avg'; % {'Asc/Desc'}
coThres         = 0.25; % coherence cutoff - between [0 1]

includeRuns     = 2;

%% setup directories
projectDir      = fullfile('/Volumes', 'server', 'Projects', projectName);
designFolder    = fullfile(projectDir,'BehavioralData');


%% Find participants
subjList        = dir(fullfile(designFolder, 'TactileSweep*'));
numSubjects     = numel(subjList);
subjNames       = {subjList.name}';

%%
for sub = 1:numSubjects
    
    strSplit        = split(subjNames{sub}, {'_' '.'});
    subject         = strSplit{contains(strSplit, 'sub')};
    
    sessions        = strSplit{contains(strSplit, 'ses-')};

    %% Directories
    dataDir         = fullfile(projectDir, 'derivatives', dataFolder,...
                        subject, sessions, 'func');
    mainDir         = fullfile(projectDir,'Data', 'phaseEncodedAnal', subject);
    saveDir         = fullfile(mainDir, 'TSeries');
    outDir          = fullfile(saveDir, 'volOverlays');
    
    if ~exist(outDir, 'dir'); mkdir(outDir); end
    
    %% See if FFT needs to be recomputed
    if isempty(dir(fullfile(outDir, 'TSeries_avg', '*.nii.gz'))) || recompute == true
        
        %% Verify fMRI data
        [session, tasks, runnums]   = bidsSpecifyEPIs(projectDir, subject(5:end),...
            sessions(5:end), taskNames(1), []);
 
        %% Get fMRI data
        if exist(fullfile(mainDir, 'Data', sprintf('%s_ses-%s_rawData.mat',subject, session)), 'file') == 0
            
            %% Setup name for JSON file containing info on the data
            fname           = sprintf('%s%s_sub-%s_ses-%s_inputVar.json', ...
                allModalities{whichModality}, allDesigns{1}, subject, session);
            inputVar        = struct('projectDir', projectDir, 'projectName', projectName, 'subject', subject, 'Modality',  allModalities{whichModality}, ...
                'Design', allDesigns{whichDesign}, 'session', session, 'tasks', taskNames, 'runnums', includeRuns, ...
                'dataFolder', dataDir, 'dataStr', dataStr, 'designFolder', designFolder, ...
                'stimdur', stimdur, 'nCycles', nCycles, 'durationCycle', durCycle, 'durationScan', durScan, 'dropFrames', dropFrames, ...
                'tr', tr, 'BOLDshift', shift, 'coherenceThreshold', coThres);
            
            savejson('',inputVar,fullfile(saveDir,fname));
            
            
            fprintf('Loading data')
            
            if ~exist(fullfile(mainDir, 'Data'), 'dir'); mkdir(fullfile(mainDir, 'Data')); end
            
            data                        = cell(size(runnums{1}));
            info                        = cell(size(runnums{1}));
            
            for jj = 1:length(runnums{1})
                fprintf('.')
                
                % we want to check for both 0-padded and non-0-padded versions...
                fnamePrefix         = sprintf('*_task-%s*run-%d_*%s*',...
                    tasks{1},runnums{1}(jj), dataStr);
                fnamePrefixZeroPad  = sprintf('*_task-%s*run-%02d_*%s*',...
                    tasks{1},runnums{1}(jj), dataStr);
                
                fname = dir(fullfile(dataDir, fnamePrefix));
                % we only need to check both if they're different; if we're
                % looking at run 10, 0-padded and non-0-padded will be the
                % same string
                if ~strcmp(fnamePrefix, fnamePrefixZeroPad)
                    fname = [fname; dir(fullfile(dataDir, fnamePrefixZeroPad))];
                end
                
                assert(length(fname) == 1);
                
                fullDatafile  = MRIread(fullfile(dataDir, fname.name));
                data{jj}      = single(fullDatafile.vol); % MRIread reorders dims
                info{jj}      = fullDatafile;
                info{jj}      = rmfield(info{jj}, 'vol');
                %{
                fullDatafile  = niftiread(fullfile(dataDir, fname.name));
                data{jj}      = single(fullDatafile);
                info{jj}      = niftiinfo(fullfile (dataDir, fname.name));
                %}
            end
            
            save(fullfile(mainDir, 'Data', sprintf('%s_ses-%s_rawData.mat',subject, session)), 'data', 'info', '-v7.3')
        else
            load(fullfile(mainDir, 'Data', sprintf('%s_ses-%s_rawData.mat',subject, session)), 'data', 'info');
            
        end
        
        %% Get behavioral data
        switch subject
            case 'sub-wlsubj114'
                upSweepRuns     = 1;
                downSweepRuns   = 2;
            case 'sub-wlsubj127'
                upSweepRuns     = 2;
                downSweepRuns   = 1;
            otherwise % When used as a localizer - only ascending sweeps
                upSweepRuns     = 1:includeRuns;
                downSweepRuns   = 0;
        end
        
        %% traveling wave analysis
        fprintf('\n Analyzing data .. \n')
        % read in volumes
        fspth               = fullfile(projectDir, 'derivatives', 'freesurfer', sprintf('sub-%s',subject));
        
        % setup boxcar kernel for smoothing / highpass filter
        smoothperiod    = 25;
        kernel          = ones([smoothperiod 1]) / smoothperiod;
        numIterations   = 2;
        
        % loop over tasks
        for kk = 1:length(tasks)
            
            % indicate task in console
            disp(tasks{kk});
            
            % read out dimensions of the data
            nFrames = info{1}.nframes;
            volDims = info{1}.volsize;
            nVoxels = prod(volDims);
            
            % initialize matrix to hold prepared, non-averaged time series
            TSeries         = NaN(numel(runnums{kk}), nVoxels, nFrames-(dropFrames*2));
            flippedTSeries  = NaN(numel(runnums{kk}), nVoxels, nFrames-(dropFrames*2));
            % setup linear model / contrast
            model           = [linspace(0,1,nFrames)' ones(nFrames,1)];
            
            
            % iterate over the runs
            for ii = runnums{kk}
                
                % reshape 3d volume into long vector: voxels x time
                reshapeData = reshape(data{ii}, [nVoxels, nFrames]);
                
                % calculate mean per voxel
                dc          = nanmean(reshapeData,2);
                % prevents divide-by-zero
                dc(dc==0 | isnan(dc)) = Inf;
                % divide each vertex's time series by its mean
                detrendData = bsxfun(@rdivide, reshapeData, dc);
                
                % Percent signal change
                detrendData = (detrendData - 1) * 100;
                
                % remove linear trend
                wgts    = model\detrendData';
                trends  = model*wgts;
                detrendData = detrendData - trends';
                
                % Initialize the baseline array with 1-period
                %{
        % padding at beginning and end:
        bLine           = zeros(nVoxels, nFrames + 2*smoothperiod);
        % set front padding to mean of period length data at the beginning of run
        bLine(:, 1:smoothperiod) = repmat(mean(reshapeData(:,1:smoothperiod), 2),1,smoothperiod);
        % set main part to run data
        bLine(:, smoothperiod+1:smoothperiod+nFrames) = reshapeData;
        % set front padding to mean of period length data at the beginning of run
        bLine(:, smoothperiod+nFrames+1:nFrames + 2*smoothperiod) = repmat(mean(reshapeData(:,1:smoothperiod), 2),1,smoothperiod);
        % Smoothing loop -- convolve with kernel for number of iterations
        for jj=1:numIterations, bLine = conv2(bLine', kernel)'; end
        % trim extra point from convolution and padding off
        bLine = bLine(:, floor(numIterations*(smoothperiod - 1)/2) + smoothperiod + 1 : floor(numIterations * (smoothperiod - 1)/2) + smoothperiod + nFrames);
        % Remove baseline from time series:
        reshapeData = reshapeData - bLine;
        % clean house
        clear('bLine');
                %}
                
                % Crop first and last 5 TRs
                TSeries(ii,:,:) = detrendData(:,dropFrames+1:nFrames-dropFrames);
                if any(ismember(upSweepRuns, ii))
                    flippedTSeries(ii,:,:) = detrendData(:,dropFrames+1:nFrames-dropFrames);
                elseif any(ismember(downSweepRuns, ii))
                    flippedTSeries(ii,:,:) = fliplr(detrendData(:,dropFrames+1+shift:nFrames-dropFrames+shift));
                end
                
            end
            
            % Average over runs
            switch whichAnalysis
                case 'Avg'
                    TSeries_avg     = squeeze(mean(flippedTSeries,1));
                    allData         = {TSeries_avg};
                    dataNames       = {'TSeries_avg'};
                case 'Asc/Desc'
                    TSeries_avg     = squeeze(mean(flippedTSeries,1));
                    TSeries_up      = squeeze(mean(TSeries(upSweepRuns,:,:),1));
                    TSeries_down    = squeeze(mean(TSeries(downSweepRuns,:,:),1));
                    allData         = {TSeries_avg, TSeries_up, TSeries_down};
                    dataNames       = {'TSeries_avg', 'TSeries_asc', 'TSeries_desc'};
            end
            
            
            for nn = 1:numel(allData)
                
                analDir = fullfile(outDir, dataNames{nn});
                if ~exist(analDir, 'dir'); mkdir(analDir); end
                
                % Compute FFT for each voxel
                %
                ft = fft(allData{nn}');
                ft = ft(1:1+fix(size(ft, 1)/2), :);
                
                % This quantity is proportional to the amplitude
                %
                scaledAmp = abs(ft);
                
                % This is in fact, the correct amplitude
                %
                amp = 2*(scaledAmp(nCycles+1,:))/size(allData{nn},2);
                
                % Compute power over the whole spectrum and calculate coherence
                %
                sqrtsummagsq = sqrt(sum(scaledAmp.^2));
                co = scaledAmp(nCycles+1,:)./sqrtsummagsq;
                
                % Calculate phase:
                % 1) add pi/2 so that it is in sine phase.
                % 2) minus sign because sin(x-phi) is shifted to the right by phi.
                % 3) Add 2pi to any negative values so phases increase from 0 to 2pi.
                %
                ph = angle(ft(nCycles+1,:));
                ph(ph<0) = ph(ph<0)+pi*2;
                
                
                % reshape parameters
                volAmp      = single(reshape(amp, [volDims 1]));
                volCo       = single(reshape(co, [volDims 1]));
                volPh       = single(reshape(ph+1, [volDims 1]));
                
                % write volumes
                saveInfo =  info{1};
                saveInfo.Description = 'Modified using MATLAB R2020b'; 
                saveInfo.nframes = 1;
                saveInfo.vol = volAmp;
                MRIwrite(saveInfo, fullfile(analDir, ['Amp_' subject '_allRuns_wholeBrain.nii.gz' ]));
                saveInfo.vol = volCo;
                MRIwrite(saveInfo, fullfile(analDir, ['Co_' subject '_allRuns_wholeBrain.nii.gz' ]));
                saveInfo.vol = volPh;
                MRIwrite(saveInfo, fullfile(analDir, ['Phase_' subject '_allRuns_wholeBrain.nii.gz' ]));

                % Thresholded volumes
                thresIndx               = volCo <= coThres;
                volThresAmp             = volAmp;
                volThresAmp(thresIndx)  = NaN;
                volThresCo              = volCo;
                volThresCo(thresIndx)   = NaN;
                volThresPh              = volPh;
                volThresPh(thresIndx)   = NaN;
                
                
                % write thresholded volumes
                saveInfo.vol = volThresAmp;
                MRIwrite(saveInfo, fullfile(analDir, ['thresAmp_' subject '_allRuns_wholeBrain.nii.gz' ]));
                saveInfo.vol = volThresCo;
                MRIwrite(saveInfo, fullfile(analDir, ['thresCo_' subject '_allRuns_wholeBrain.nii.gz' ]));
                saveInfo.vol = volThresPh;
                MRIwrite(saveInfo, fullfile(analDir, ['thresPhase_' subject '_allRuns_wholeBrain.nii.gz' ]));

                
                % Write averaged timeserie volume
                volTSeries = reshape(TSeries_avg, [volDims size(TSeries_avg,2)]);

                saveInfo.nframes = size(TSeries_avg,2);
                saveInfo.vol = volTSeries;
                MRIwrite(saveInfo, fullfile(mainDir, 'Data', sprintf('sub-%s_ses-%s_avgTSeries_%.0fRuns.nii.gz', subject, session, includeRuns)));
                
                
            end
        end
    end
end
