% function to re-format the results multcompare, using the stats from the 
% one-way anova function, into a more readable format

function mc = multcompare_ph(stats)
    
    % extract the group names from stats
    group_names = stats.gnames;
    
    % perform post-hoc comparison
    [results,means] = multcompare(stats);
    
    % add in the group names to the mc results, if specified
    if (~isempty(group_names))
        column1_gnames = cell(size(results,1),1);
        column2_gnames = cell(size(results,1),1);
        for i=1:size(results,1)
            column1_gnames{i} = group_names{results(i,1)};
            column2_gnames{i} = group_names{results(i,2)};
        end     
    else
        column1_gnames = results(:,1);
        column2_gnames = results(:,2);
    end

    mc = table(column1_gnames, column2_gnames, results(:,3),results(:,4),results(:,5),results(:,6),...
                'VariableNames', {'group1', 'group2', 'g1_minus_g2_95CI_low', 'g1_minus_g2_mean', 'g1_minus_g2_95CI_upp', 'pvalue'});

end