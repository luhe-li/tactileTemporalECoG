% Apply Bipolar/Laplacian referencing per shaft for specified subjects and tasks
rootPath = '/Volumes/Epilepsy_ECOG/SharedAnalysis';
addpath('/Volumes/Epilepsy_ECOG/CodeBase/Matlab/');
close all;

% --- Define Subjects and Tasks to Process ---
SUBJs = {'NY905'};
tasks = {'VisRead'};

% --- Main Loop over Subjects ---
for j = 1:numel(SUBJs)
    SJ = SUBJs{j};
    fprintf('Current subject is %s\n', SJ);
    
    % Get globals to find coordinates file
    get_subj_globals(SJ,'PicN',rootPath); 
    
    % Read coordinates to define the start and end of each shaft
    T = readtable([SJdir dlm '/analysis/coordinates.csv']);
    labelLetter=string;
    for i=1:height(T)
        current_loc_pos = regexp(T.labels(i,1),'\d*','split');
        labelLetter(i,1) = current_loc_pos{1}{1,1};
    end
    
    shaftEnd = [find(~strcmp(labelLetter(1:end-1), labelLetter(2:end))); length(labelLetter)];
    shaftStart = [1; shaftEnd(1:end-1)+1];
    
    % Find electrodes that don't exist (e.g., NaN coordinates)
    elecs_nonexist = find(isnan(T.T1_x));
     fprintf('There are %d electrodes that are outside of the shaft \n', numel(elecs_nonexist));
    
    % --- Loop over Tasks for the Current Subject ---
    for t = 1:numel(tasks)
        task = tasks{t};
        fprintf('Processing Task: %s\n', task);
        get_subj_globals(SJ,task,rootPath);
        
        % Load raw data (gdat) and bad electrodes (bad_elecs)
        load([DTdir dlm 'gdat.mat']) 
        
        % Demean all electrodes first
        gdat_demean = gdat - mean(gdat,2);
        
        % Initialize the output matrix
        gdat_referenced = gdat_demean; 
        
        % --- Loop over Electrode Shafts ---
        for i = 1:length(shaftStart)
            elec_crt = shaftStart(i):shaftEnd(i); % Electrodes in the current shaft
            elec_car_crt = elec_crt; % Copy for modification
            
            % Remove bad electrodes from the referencing pool
            if sum(ismember(elec_crt,bad_elecs)~=0)
                elec_car_crt(find(ismember(elec_car_crt,bad_elecs)~=0)) = [];
            end
            
            % Remove non-existing electrodes from the referencing pool
            if sum(ismember(elec_crt,elecs_nonexist) ~=0)
                elec_car_crt(find(ismember(elec_car_crt,elecs_nonexist)~=0)) = [];
            end
        
            % Only proceed if there are good electrodes left for referencing
            if ~isempty(elec_car_crt)
                gdat_crt = gdat_demean(elec_car_crt,:);
                gdat_CAR_crt = zeros(size(gdat_crt)); % Temp matrix for referenced data
                
                shaft_num_elec = size(gdat_CAR_crt,1);

                if shaft_num_elec == 1
                    % If only one good electrode, it keeps its demeaned value
                    gdat_CAR_crt(1,:) = gdat_crt(1,:);
                
                elseif shaft_num_elec > 1
                    % First elec is referenced bipolarly to the next one
                    gdat_CAR_crt(1,:) = gdat_crt(1,:) - gdat_crt(2,:);
                
                    % Last elec is referenced bipolarly to the previous one
                    gdat_CAR_crt(shaft_num_elec,:) = gdat_crt(shaft_num_elec,:) - gdat_crt(shaft_num_elec-1,:);
                
                    % Inner electrodes are referenced using Laplacian
                    % (2*current - (prev + next))
                    for k = 2:shaft_num_elec-1
                        gdat_CAR_crt(k,:) = gdat_crt(k,:).*2 - (gdat_crt(k-1,:) + gdat_crt(k+1,:));
                    end
                end
                
                % Put the referenced data back into the main matrix
                for k = 1:numel(elec_car_crt)
                    gdat_referenced(elec_car_crt(k),:) = gdat_CAR_crt(k,:);
                end
                
                fprintf('Finished shaft %d. (%d good elecs)\n',i, shaft_num_elec)
            
            else % No good electrodes in the shaft
                % Keep the original demeaned data for these electrodes
                gdat_referenced(elec_crt,:) = gdat_demean(elec_crt,:);
                fprintf('Shaft %d has no good electrodes. Using demeaned data.\n',i)
            end
        end
        
        % Overwrite gdat with the newly referenced data
        gdat = gdat_referenced;
        
        % Save the referenced data
        save([DTdir dlm 'gdat_bipolar.mat'],"gdat")
        fprintf('Saving gdat_bipolar.mat for task %s...\n', task)
    end % end task loop
end % end subject loop

fprintf('All subjects and tasks processed.\n');