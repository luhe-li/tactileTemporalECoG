function tt_writeParamsToJson(modelName)
% Writes out a json file with parameters for model fitting.
% 2020 Iris Groen

%% define starting points and bounds on parameters

switch modelName
    
    case {'DN', 'DNCASCADE'}
        
        opts.params = 't1,w,t2,n,sigma,shift,scale';
        opts.x0     = [0.03, 0,   0.07, 1.5, 0.15, 0.06, 2];    % starting point
        opts.lb     = [0.01, 0,   0.01, 1,   0,    0,    0.01]; % lower bounds
        opts.ub     = [1,    1,   2,    5,   1,    0.1,  200];  % upper bounds
        opts.plb    = [0.1,  0,   0.1,  1.5, 0.01, 0.01, 0.5];  % plausible lower bound (required for bads search algorithm)
        opts.pub    = [0.9,  0.5, 1,    3,   0.5,  0.08, 100];  % plausible upper bound (required for bads search algorithm)
    
	case {'LINEAR', 'LINEAR_RECTH', 'LINEAR_RECTF'} 
         
        opts.params = 'tau1, tau2, n_irf, weight, shift, scale';
        opts.x0     = [0.01, 0.03, 2,     0,      0.01,  2];    % starting point
        opts.lb     = [0.01, 0.01, 2,     0,      0,     0.01]; % lower bounds
        opts.ub     = [10,   10,   20,    1,      1,     200];  % upper bounds
        opts.plb    = [0.1,  0.1,  2,     0.01,   0.01,  0.5];  % plausible lower bound (required for bads search algorithm)
        opts.pub    = [0.9,  0.9,  10,    0.5,    0.08,  100];  % plausible upper bound (required for bads search algorithm)
    
end

%% write out json
fname = fullfile(tt_RootPath, 'models', sprintf('%s.json', modelName));
savejson('',opts,fname);