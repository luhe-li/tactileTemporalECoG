function [stim, stim_info] = tt_simulateNewStimuli(t,nStim)

% Generate stimulus time courses to use for generating new model
% predictions after having fitted a model to the experimental data.
%
% 2022 Iris Groen

% DUR
% generate stimuli with nStim durations between 0 and 1.2;
stimDUR = zeros(length(t),nStim);
nameDUR = cell(nStim,1);
durationDUR = zeros(nStim,1);
isiDUR = zeros(nStim,1);

for ii = 1:nStim
    t_off = ii/nStim * 1.2;
    stim_on = t>0 & t<=t_off;
    stimDUR(stim_on,ii) = 1;
    durationDUR(ii) = round(t_off,3);
    nameDUR{ii} = sprintf('ONEPULSE-sim%d', ii);
end

% ISI
% generate stimuli with nStim ISIs between 0 and 0.533;
stimISI = zeros(length(t),nStim);
nameISI = cell(nStim,1);
durationISI = ones(nStim,1)*0.4;
isiISI = zeros(nStim,1);
% generate first pulse
pulse1_toff = 0.2;
pulse1_on = t>0 & t<=pulse1_toff;
stimISI(pulse1_on,:) = 1;
ISIs = linspace(0,1.2,nStim);
% generate second pulse
for ii = 1:nStim
    this_ISI = ISIs(ii);
    pulse2_ton = pulse1_toff + this_ISI;
    pulse2_toff = pulse2_ton + 0.2;
    pulse2_on = t>pulse2_ton & t<=pulse2_toff;
    stimISI(pulse2_on,ii) = 1;
    isiISI(ii) = round(this_ISI,3);
    nameISI{ii} = sprintf('TWOPULSE-sim%d', ii);
end

% concatenate conditions
stim = horzcat(stimDUR, stimISI);
name = vertcat(nameDUR, nameISI);
duration = vertcat(durationDUR, durationISI);
ISI =  vertcat(isiDUR, isiISI);

% generate new stim_info table
stim_info = table(name, duration, ISI);    

end