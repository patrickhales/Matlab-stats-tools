% ttest_ph.m (c) Patrick Hales, 2017. 
% Tested on:
% Matlab R2016a
% NOT compatible with:
% Matlab R2013a
%
% Function to compare two samples (x and y). This will automatically check the
% normality of the samples, and perform the appropriate parametric / non-parametric test. 
% The null hypothesis is that x and y are samples drawn from a polulation with the same
% mean / median (depending on the test used). 
%
% Tests used are...
% -- Parametric tests (normally distributed data)
%       - Un-paired:    Two-tailed, two-sample t-test (Behrens-Fisher correction for unequal variances will be used if needed)
%       - Paired:       Two-tailed paired t-test
% -- Non-parametric tests (not normally distributed data)
%       - Un-paired:    Two-sided Wilcoxon rank sum test (same as Mann-Whitney U-test)
%       - Paired:       Wilcoxon signed rank test
%
% Required arguments:
%   x = sample 1 (1D vector of data)
%   y = sample 2 (1D vector of data)
% Output:
%   h = 0 : failure to reject the null hypothesis (at the significance level - default is 5%)
%           (i.e. samples drawn from populations with the same mean/median)
%   h = 1 : rejection of the null hypothesis (at the significance level)
%           (i.e. samples drawn from populations with different mean/median)
%   p =     p-value 
%
% Examples of use:
% 
% (1) default settings: performs un-paired test at default 5% significance level
% [h, p, test_flag] = ttest_ph(x,y)  
%
% (2) perform a paired test (note x and y must be the same length)
% [h, p, test_flag] = ttest_ph(x,y,'paired')
% (when needed, 'paired' must always be the 3rd argument. The rest of the 
% name/value pair arguments below can appear anywhere)
%
% (3) specify a significance level of, say, 1%
% [h, p, test_flag] = ttest_ph(x,y,'alpha',0.01)
%
% (4) Show histograms and boxplots of the two samples
% [h, p, test_flag] = ttest_ph(x,y,'plot', 'on')
%
% Notes on the tests used:
%
% Using Lilliefors test to check normality:
% Null hypothesis is data in vector x comes from a distribution in the 
% normal family
% The result h is 1 if the test rejects the null hypothesis at the 5% 
% significance level, and 0 otherwise.
% So    h=0: data is normal
%       h=1: data is not normal
% If you lower the significance level, it is harder to reject the null
% hypothesis that the data is normal (using e.g. 'Alpha',0.01)
%
% Using F-test to check for equal variance of two groups (Two-sample F-test):
% Null hypothesis is that the data in vectors x and y comes from normal 
% distributions with the same variance. The alternative hypothesis is 
% that they come from normal distributions with different variances.
% The result h is 1 if the test rejects the null hypothesis at the 5% 
% significance level, so..
%       h=0: x and y have equal variances
%       h=1: x and y have different variances


function [h, p, test_flag] = ttest_ph(x,y,varargin)

    % Default settings:
    defaultPaired = 'unpaired';
    defaultAlpha = 0.05;
    defaultPlots = 'off';
    defaultXcolour = 'b';
    defaultYcolour = 'g';
    
    args = inputParser;
    addRequired(args,'x',@isnumeric);
    addRequired(args,'y',@isnumeric);   
    validFlag1 = {'paired'};
    checkFlag1 = @(x) any(validatestring(x,validFlag1));
    addOptional(args, 'paired_flag', defaultPaired, checkFlag1);
    addParameter(args, 'alpha_user', defaultAlpha, @isnumeric);
    validOnOff = {'on','off'};
    checkValidOnOff = @(x) any(validatestring(x,validOnOff));
    addParameter(args, 'plots', defaultPlots, checkValidOnOff);
    % each parameter added with addParameter corresponds to two input 
    % arguments: one for the name and one for the value of the parameter.
    % addOPtional only required one argument
    parse(args,x,y,varargin{:});

    % check input data
    if (~isvector(x) || ~isvector(y))
        error('Pass two 1D vectors');
    end
    % make each vector a row vector
    if (size(x,1) > 1)
        x = x';
    end
    if (size(y,1) > 1)
        y = y';
    end
    
    if (strcmp(args.Results.paired_flag, 'paired'))
        paired_flag = 1;
    else
        paired_flag = 0;
    end
    alpha_user = args.Results.alpha_user;
    if (paired_flag && (length(x) ~= length(y)))
        error('x and y must be equal lengths for paired test');
    end
    plots_flag = args.Results.plots;
       
    disp(' ');
    if (paired_flag==0)
        % ---- Un-paired tests --------------------------------------------
        % For un-paired tests, we need to check data in each group is
        % normally distributed
        [x_not_norm, px] = lillietest(x);
        [y_not_norm, py] = lillietest(y);
        disp(' ');
        % we also need to check if the two groups have equal variance
        [non_equal_var, pvar] = vartest2(x,y);
    
        if (x_not_norm)
            disp(['   x is NOT normal (p=' num2str(px) ' for rejecting normality)']);
        else
            disp(['   x is normal (p=' num2str(px) ' for rejecting normality)']);
        end
        if (y_not_norm)
            disp(['   y is NOT normal (p=' num2str(py) ' for rejecting normality)']);
        else
            disp(['   y is normal (p=' num2str(py) ' for rejecting normality)']);
        end
            
        if ((x_not_norm==0) & (y_not_norm==0))
            % Parametric tests
            if (non_equal_var==0)
                disp('Performing Students t test with equal variances...');
                test_flag = 'TtestEqualVar';
                [h,p] = ttest2(x,y,'Alpha',alpha_user);
            else
                disp('Performing Students t test with non-equal variances (Behrens-Fisher)...');
                test_flag = 'TtestNonequalVar';
                [h,p] = ttest2(x,y,'Vartype','unequal','Alpha',alpha_user);
            end
        else   
            % Non-parametric tests
            disp('Performing two-sided Wilcoxon rank sum test (Mann-Whitney U-test)...')
            test_flag = 'WilcoxonRankSum';
            [p,h] = ranksum(x,y,'Alpha',alpha_user);
        end
        % -----------------------------------------------------------------
    else
        % -------- Paired tests -------------------------
        % For paired tests, we are checking that 'x-y' is normally
        % distributed (and possibly not centred on zero). Equal variance is
        % not an issue here, as we are only looking at one distribution
        % (x-y).
        groups_diff = x - y;
        [gdiff_not_norm, pgdiff] = lillietest(groups_diff);
                
        if (gdiff_not_norm==0)
            % Parametric tests
            disp(['Difference between x and y is normally distributed (p=' num2str(pgdiff) ')']);
            disp('Performing paired Students t test...');
            test_flag = 'TtestPaired';
            [h,p] = ttest(x,y,'Alpha',alpha_user);
        else   
            % Non-parametric tests
            disp(['Difference between x and y is NOT normally distributed (p=' num2str(pgdiff) ')']);
            disp('Performing Wilcoxon signed rank test...')
            test_flag = 'WilcoxonSignedRank';
            [p,h] = signrank(x,y,'Alpha',alpha_user);
        end
        % -----------------------------------------------------------------
    end
        
    % General
    disp(' ');
    disp(['-> Median X = ' num2str(median(x(~isnan(x))))]);
    disp(['-> Median Y = ' num2str(median(y(~isnan(y))))]);
    if (h==1)
        disp(['-> Significant difference (at alpha = ' num2str(alpha_user) ') between x and y (p=' num2str(p) ')']);
    else
        disp(['-> No difference (at alpha = ' num2str(alpha_user) ') between x and y (p=' num2str(p) ')']);
    end
    disp(' ');
    
    % Plot results (if requested)
    if (strcmp(plots_flag, 'on'))
       close all; % tmp!
       figure()
       subplot(1,2,1);
       if (length(x) > 10 && length(y) > 10)
           h1 = histogram(x,'Normalization','probability', 'FaceColor', defaultXcolour, 'FaceAlpha', 0.5);
           hold on
           h2 = histogram(y,'Normalization','probability', 'FaceColor', defaultYcolour, 'FaceAlpha', 0.5);
       else
           h1 = histogram(x,'Normalization','probability', 'FaceColor', defaultXcolour, 'FaceAlpha', 0.5, 'NumBins', 10);
           hold on
           h2 = histogram(y,'Normalization','probability', 'FaceColor', defaultYcolour, 'FaceAlpha', 0.5, 'NumBins', 10);
       end           
              
       if (h1.NumBins > h2.NumBins)
           h2.NumBins = h1.NumBins;
       end
       if (h2.NumBins > h1.NumBins)
           h1.NumBins = h2.NumBins;
       end
       if (h1.BinWidth < h2.BinWidth)
           h2.BinWidth = h1.BinWidth;
       end
       if (h2.BinWidth < h1.BinWidth)
           h1.BinWidth = h2.BinWidth;
       end
   
       subplot(1,2,2);
       jitterplot_ph_ttest_ph([x y],[zeros(1,length(x)) ones(1,length(y))]);
           
    end  
end

% -------------------------------------------------------------------------
% Jitterplot function for plotting:
function jitterplot_ph_ttest_ph(y,x,jitter)

    figure(gcf)
    lw = 1.5;   % default line width
    fs = 15;    % default font size
    
    
    if (isvector(y))
        y=y(:); 
    end
    if (nargin<2 || isempty(x))
        x=1:size(y,2);
    end
    
    if (isvector(y) && isvector(x) && length(x)>1)
        x=x(:);
        if length(x)~=length(y)
            error('length(x) should equal length(y)')
        end
        u=unique(x);    % u is a list of group numbers
    end
    
    if nargin<3 || isempty(jitter)
        jitter=0.3;
    end
    h0 = boxplot(y, x, 'colors', 'k', 'symbol', 'k');
    hold on
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),[0.85 0.85 0.85]);
        set(h0(j,:),'LineWidth',lw);
    end
    set(gca,'XTickLabel',{' '})
    delete(h0);
    h0 = boxplot(y, x, 'colors', 'k', 'symbol', 'b');
    for j=1:length(h0)
        set(h0(j,:),'LineWidth',lw);
    end
    set(gca,'XTickLabel',{' '})
    
    xvals_all0 = [];
    yvals_all0 = [];
    rawinds_list = [];
    
    % loop through groups
    for i=1:length(u)
        if (i==1)
            figure(gcf);
            hold on
        end
        f=find(x==u(i));
        ygroup = y(f);
        xbase = i; %u(i);
        rawinds_list = [rawinds_list f'];
        xgroup_jittered = xbase + (rand(size(ygroup))-0.5)*jitter;
        ygroup_mean = mean(ygroup);
        ygroup_std = std(ygroup);
        minx = xbase - min(xgroup_jittered);
        maxx = max(xgroup_jittered) - xbase;
        reach=0.10;

        hold on
        plot(xgroup_jittered, ygroup, 'ko', 'markerfacecolor', [0.8 0.8 1.0]);
        
        xvals_all0 = [xvals_all0 xgroup_jittered'];
        yvals_all0 = [yvals_all0 ygroup'];

    end

    set(gca, 'linewidth',lw, 'FontSize',fs)
    hold off
    [vals, inds] = sort(rawinds_list);
    xvals_all = xvals_all0(inds);
    yvals_all = yvals_all0(inds);
end
% -------------------------------------------------------------------------
     