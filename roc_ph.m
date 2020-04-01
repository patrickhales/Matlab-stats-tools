% function to calculate ROC curves. Note this is built into Matlab, but it's confusing!
% Pass in two vectors of values (one for each group). The optimum threshold
% for splitting these (based on best separation of the two groups) will be found

% required arguments: 
%   y1, y2: two vectors of data, one from each group

% optional arguments:
%   'plot', 'on':       show plots of ROC curve and jitterplot of two
%   groups

% simulated data for testing:
    % npoints=40;
    % y1 = normrnd(20,2,[1,npoints/2]);
    % y2 = normrnd(10,2,[1,npoints/2]);


function [optimum_threshold, max_sensitity, max_specificity, area_under_curve] = roc_ph(y1, y2, varargin)
    
    % default settings
    defaultPlots = 'off';
    
    args = inputParser;
    addRequired(args,'y1',@isnumeric);
    addRequired(args,'y2',@isnumeric);  
    
    validOnOff = {'on','off'};
    checkValidOnOff = @(x) any(validatestring(x,validOnOff));
    addParameter(args, 'plots', defaultPlots, checkValidOnOff);
    parse(args,y1,y2,varargin{:});
    plots_flag = args.Results.plots;

    % check data are row vectors
    if (ndims(y1) > 2) | (ndims(y2) > 2)
        error('Pass two 1D data vectors');
    end
    if (size(y1,1) > 1)
        y1 = y1';
    end
    if (size(y2,1) > 1)
        y2 = y2';
    end
    
    groups = [zeros(1,numel(y1)), ones(1,numel(y2))];
    y = [y1 y2];

    total_actual_positives = sum(groups==1);
    total_actual_negatives = sum(groups==0);

    steps = 1000;
    n = numel(y);

    splitvec = linspace(min(y), max(y), steps); 

    mean_group0 = mean(y(groups==0));
    mean_group1 = mean(y(groups==1));

    if (mean_group0 >= mean_group1)
        g0high = 1;
    else
        g0high=0;
    end

    predicted_groups = zeros(steps,n);
    true_positive_rate = zeros(1,steps);
    false_positive_rate = zeros(1,steps);

    for t=1:steps
        thresh = splitvec(t);
        if (g0high)
            predicted_groups(t,y<thresh) = 1;
        else
            predicted_groups(t,y>thresh) = 1;
        end
        truepositives_vec = (squeeze(predicted_groups(t,:)) == 1) & (groups==1);
        truenegatives_vec = (squeeze(predicted_groups(t,:)) == 0) & (groups==0);
        falsepositives_vec = (squeeze(predicted_groups(t,:)) == 1) & (groups==0);
        falsenegatives_vec = (squeeze(predicted_groups(t,:)) == 0) & (groups==1);

        true_positive_rate(t) = sum(truepositives_vec) / total_actual_positives;    % sensitivity
        false_positive_rate(t) = sum(falsepositives_vec) / total_actual_negatives;  % 1 - specificity
    end
    sensitivity = true_positive_rate;
    specificity = 1 - false_positive_rate;

    distance_from_top_left = sqrt( (false_positive_rate.^2) + ((1-true_positive_rate).^2) );
    [val, optimum_ind] = min(distance_from_top_left);
    area_under_curve = abs(trapz(1-false_positive_rate, true_positive_rate));
    % [val, optimum_ind] = max(sensitivity+specificity);
    optimum_threshold = splitvec(optimum_ind);
    max_sensitity = sensitivity(optimum_ind);
    max_specificity = specificity(optimum_ind);
    
    if (strcmp(plots_flag, 'on'))
        figure()
        subplot(1,2,1);
        plot([0,1],[0,1], 'k--');
        hold on
        plot(false_positive_rate, true_positive_rate, 'k-');
        xlabel('False positive rate');
        ylabel('True positive rate');
        subplot(1,2,2)
        jitterplot_ph(y,groups);
        hold on
        xr = xlim;
        plot([xr(1),xr(2)], [optimum_threshold, optimum_threshold], 'k--');
    end

%     disp(['Threshold = ' num2str(optimum_threshold) ':    Max sensitivity = ' num2str(max_sensitity) ', max specificity = ' num2str(max_specificity), ' AUC = ' num2str(area_under_curve)]);
%     [h,p] = ttest2(y(groups==0), y(groups==1));
%     num2clip([optimum_threshold, max_sensitity, max_specificity, p]);

end

    