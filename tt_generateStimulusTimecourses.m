function [stim, stimInfo] = tt_generateStimulusTimecourses(stimNames,t, stimInfo)

% Generates stimulus timecourses based on stimulus name and epoch duration,
% provided in a stimInfo table; will load stimInfo table in the code
% repository if no stimInfo input argument is provided.

% Rerote from tde_generateStimulusTimecourses, fixed contrast to 1

% <subjectList>
if ~exist('stimInfo', 'var') || isempty(stimInfo)
    fprintf('[%s] Loading stimulus info \n',mfilename);
    stimInfo_fname = fullfile(tt_RootPath, 'stiminfo_ny.tsv');
    stimInfo = readtable(stimInfo_fname, 'FileType', 'text');
end

fprintf('[%s] Generating stimulus timecourses \n',mfilename);

condition = stimInfo.name;
duration  = stimInfo.duration;
ISI       = stimInfo.ISI;
contrast  = 1;

nStim = length(stimNames);
stim  = zeros(length(t),nStim);

for cond = 1:length(stimNames)
    inx         = contains(condition, stimNames{cond});
    stim(t>0 & t<=duration(inx),cond) = contrast;
    if ISI(inx) > 0
        t_pulse1off = duration(inx)+ISI(inx);
        stim(t>t_pulse1off & t<=(t_pulse1off+duration(inx)),cond) = contrast;
    end
end
fprintf('[%s] Done! \n',mfilename);

end