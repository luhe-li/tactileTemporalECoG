function [results] = tt_plotTactileVisualParams(results, modelfun, saveDir)

% Plot parameter estimaets from tactile and visual data from the same model
% (thus the same parameters)

if ~exist('saveDir', 'var') || isempty(saveDir), saveDir = []; end
if ~isempty(saveDir), saveFig = true; else, saveFig = false; end

v_ind = 3;

clut = {[1, 0, 0],[0.7, 0.7, 0.7],[0, 0, 0]};
ymax = [0.1, 0.6, 0.4, 1.8, 0.1, 0.15, 4];

%% Prep

% process visual data not averaged
nModels     = size(results,2);
[~, v_channels, group_prob] = groupElecsByVisualArea(results{v_ind}.channels, 'probabilisticresample');
v_nChans = height(v_channels);

figure('Name', 'tact_vis_fitted parameters'); hold on
set(gcf, 'Position', [400 800 2000 600]);

%% Plot

% Read in parameter names from json
tmp = loadjson(fullfile(tdeRootPath, 'temporal_models', sprintf('%s.json', func2str(modelfun))));
paramNames = strsplit(tmp.params,',');
nParams = length(paramNames); % shared between tact and vis

for p = 1:nParams
    subplot(2,ceil(nParams/2),p); hold on

    for kk = nModels:-1:1

        % Visual (individualelecs)
        if kk == v_ind

            [m, se] = averageWithinArea(results{kk}.params(p,:), group_prob);
            errorbar(1:v_nChans, m, m-se(:,1)', se(:,2)'-m, '.','Color',clut{kk}, 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
            tks = v_nChans;

            % Tactile (electrodeaverages, individualelectrodes)
        else

            t_nChans = size(results{kk}.params,2);
            scatter(tks+1:tks+t_nChans, results{kk}.params(p,:), 60,  clut{kk}, 'filled')
            tks = tks + t_nChans;

            if kk == 1
                tks_label = [tks_label; {'Avg'}];
            else
                tks_label = [v_channels.name; results{kk}.channels.name];
            end

        end

    end

    set(gca, 'Xlim', [0 tks+1], 'XTick', 1:tks, 'XTickLabel', tks_label, 'XTickLabelRotation', 45);
    title(paramNames{p}); ylabel('parameter value');set(gca, 'fontsize', 16);
    set(gca, 'Ylim', [0, ymax(p)])
end

if saveFig
    if ~exist(saveDir), mkdir(saveDir); end
    figName = sprintf('tactVisFittedParams2_%s', func2str(modelfun));
    saveas(gcf, fullfile(saveDir, figName), 'png');
end

end
