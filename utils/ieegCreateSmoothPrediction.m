
function out = ieegCreateSmoothPrediction(params, model)

%-- for which model to create predictions
currModel           = str2func(sprintf('%smodel', model));

%-- create finer sampling and stimuli names 
opt                 = [];
opt.numConditions   = 482;
opt.stimdur         = linspace(0, 1.2, opt.numConditions/2);
opt.x_data          = 0:1:10;
opt.twoPulseDur     = 0.2;

onePulseNames       = strrep('ONE-PULSE-%s','%s', cellstr(string((1:opt.numConditions/2)')));
twoPulseNames       = strrep('TWO-PULSE-%s','%s', cellstr(string((1:opt.numConditions/2)')));

opt.ConditionNames  = cat(1, onePulseNames, twoPulseNames); 

%-- create prediction
pred                = currModel(params, opt); 

out.pred            = pred';
out.stimdur         = opt.stimdur;