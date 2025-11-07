
% -- Visualize model components

stimdur         = 1;

numTRs          = 11;
x_data          = 0:1:numTRs-1; 
samples         = 1000;
dt              = 1/samples; % in s

finer_t         = 0:1:4*samples; % model out to 4s (max duration is 1.6s)
t_length        = numTRs * samples;

normSum         = @(x) x./sum(x(:));

% Model parameters
tau             = 0.1;
sigma           = 0.2;
gain            = 1;
p1              = 4;
p2              = 7;
w               = 0.5;
n               = 2; % fixed to 2
tau2            = 0.2;

gray = repmat(0.5, [1,3]);

projectDir = fileparts(pwd);
figDir     = fullfile(projectDir, 'figures', 'modelSchematic');
if ~exist(figDir, 'dir'), mkdir(figDir); end

%%
fig = figure('Color', [1 1 1], 'Position', [30 300 800 200]);
set(fig,'Units', 'Pixels', 'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[800 390])

%Stimulus sequence
stimSeq         = finer_t > 0 & finer_t <= stimdur * samples;

subplot(2,6,1)
plot(finer_t, stimSeq, 'Color', gray, 'LineWidth', 1)
box off; axis off
xlim([0 2*samples])
title('Stimulus', 'FontSize', 14)

% set up h1, the impulse response function
h1_t            = x_data(1):dt:2;
h1              = gampdf(h1_t, 2, tau); % assume weight = 0;
h1              = normSum(h1);

subplot(2,6,2)
plot(h1, 'k', 'LineWidth', 2)
box off; axis off
xlim([0 2*samples])
title('h1', 'FontSize', 14)

% convolve with irf to create neural predictions
linResp         = conv(stimSeq, h1, 'full');
linResp         = linResp(1:length(finer_t));
numResp         = abs(linResp).^n;

subplot(2,6,3)
plot(finer_t, numResp, 'k', 'LineWidth', 2)
hold on,
plot(finer_t, stimSeq, 'Color', gray, 'LineWidth', 1)
box off; axis off
xlim([0 2*samples])
title('num resp^2', 'FontSize', 14)

% convolve with delayed response
h2_t = x_data(1):dt:2;
h2              = normSum(exp(-h2_t / tau2));
poolResp        = conv(numResp, h2, 'full');
poolResp        = poolResp(1:length(finer_t));
demResp         = sigma.^n + abs(poolResp).^n;

subplot(2,6,4)
plot(h2, 'k', 'LineWidth', 2)
box off; axis off
% xlim([0 2*samples])
title('h2', 'FontSize', 14)

subplot(2,6,5)
plot(finer_t, demResp, 'k', 'LineWidth', 2)
hold on,
plot(finer_t, stimSeq, 'Color', gray, 'LineWidth', 1)
box off; axis off
xlim([0 2*samples])
title('dem resp', 'FontSize', 14)

% normalization response
normResp        = numResp ./ demResp;

subplot(2,6,6)
plot(normResp, 'k', 'LineWidth', 2)
hold on,
plot(finer_t, stimSeq * max(normResp), 'Color', gray, 'LineWidth', 1)
xlim([0 2*samples])
box off; axis off
title('norm resp', 'FontSize', 14)

print(gcf, fullfile(figDir, sprintf('delayedNormalization_modelSchematic')), '-dpdf')

