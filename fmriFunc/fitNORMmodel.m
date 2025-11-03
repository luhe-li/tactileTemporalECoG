function out = fitNORMmodel(freeParams, fixedParams)

mode         = fixedParams{1};

if length(fixedParams) < 5
    fixedParams{5} = 'doubleGamma';
end
switch mode

    case 'initialize'

        % init(1):  tau1, neural IRF peak time, in second
        % init(2):  sigma, semi-saturation constant
        % init(3):  response gain of the predicted BOLD responses
        % init(4):  gamma1 time to peak, in second
        % init(5):  gamma2 time to peak, in second
        % init(6):  gamma2 gain
        switch fixedParams{5}
            case 'singleGamma'
                out.opts     = optimset('display','off');
                out.init     = [  0.1     0.1       1       5 ]; % model start params
                out.lb       = [0.001    0.01       0       0 ]; % lower bound
                out.ub       = [    1       1     Inf     Inf ]; % upper bound
                out.labels   = {'tau1','sigma','gain','gamma1'};
            case 'doubleGamma'
                out.opts     = optimset('display','off');
                out.init     = [  0.1     0.1       1       5       8       1]; % model start params
                out.lb       = [0.001    0.01       0       0       0       0]; % lower bound
                out.ub       = [    1       1     Inf      20      20     Inf]; % upper bound
                out.labels   = {'tau1','sigma','gain','gamma1', 'gamma2','gamma2_gain'};
            case 'bads'
                out.opts     = optimset('display','off');
                out.init     = [  0.1     0.1       1       5       8       1]; % model start params
                out.lb       = [0.001    0.01       0       0       0       0]; % lower bound
                out.ub       = [    1       1     1e3      20      20     1e3]; % upper bound
                out.plb      = [  0.1     0.1       1       4       5     0.8]; % plausible lower bound
                out.pub      = [  0.4     0.5      20       8      10     1.2]; % plausible upper bound
                out.labels   = {'tau1','sigma','gain','gamma1', 'gamma2','gamma2_gain'};
        end
        
    case{'optimize', 'prediction'}

        x_data       = fixedParams{2};
        y_data       = fixedParams{3};
        fixedP       = fixedParams{4};

        y_est        = NORMmodel(freeParams, fixedP);
        
        assert(isequal(size(y_data,1), size(y_est,1)), 'dimensions of data and prediction are not the same')

        if strcmp(mode, 'optimize')

            if any(isnan(y_est))
                out             = 1e10;
            else
                out             = double(sum((y_data(:) - y_est(:)).^2));
            end

        elseif strcmp(mode, 'prediction')
            if ~isempty(y_data)
                out.R2      = 1 - sum((y_data(:) - y_est(:)).^2, 1) ./ sum((y_data(:) - mean(y_data(:), 1)).^2, 1);
                out.SSE     = sum((y_data(:) - y_est(:)).^2);
            end
            out.param       = freeParams;
            out.tau2        = 0;
            out.y_est       = y_est;
            out.y_data      = y_data;
            out.x_data      = x_data;
            out.condNames   = fixedP.ConditionNames;
            out.labels      = {'tau1', 'sigma','gain','gamma1', 'gamma2','gamma2_gain'};

        end

end


end