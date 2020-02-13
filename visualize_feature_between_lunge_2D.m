load('flymat_info_struct');
exp_of_interest = [3, 4];
genotype_of_interest = {'WT_SM', 'CG3385_KO_CS-bc-11_GM'};
feature_of_interest = {'facing_angle_in_lunge_interval'};
stat_of_interest = {'timepoint'};
period_of_interest = [1, 30];
is_attacked_fly_of_interest = false; 
marker_opts = {'o', 'x'};
fps = 60; 
length_of_exp = 30; 

for i=1:length(feature_of_interest)
    figure();
    hold on; 
    
    legends = {}; 
    for j=1:length(exp_of_interest)
        exp_info = flymat_info_struct(exp_of_interest(j));
        genotype = genotype_of_interest{j};
        
        if exp_info.is_same_genotype_pairing
            same_genotype_str = 'same_genotype';
        else
            same_genotype_str = 'diff_genotype';
        end
        load(strcat('featAll_and_feat_probs_dist_remove_outliers-', ...
            same_genotype_str, '-', strjoin(exp_info.flymat, '_'), '.mat'), ...
            'featAll_cell', 'lunge_interval_all_cell', ...
            'feature_options', 'stat_options', 'period_options', 'is_attacked_fly', ...
            'lunge_start_all_cell', 'num_of_flies');
        fprintf('%d %d\n', num_of_flies(1), num_of_flies(2));

        feat_mask = strcmp(feature_options, feature_of_interest{i});
        stat_mask = strcmp(stat_options, stat_of_interest{i});
        period_mask = cellfun(@(period_pair) isequal(period_pair, period_of_interest), period_options); 
        is_attacked_mask = is_attacked_fly == is_attacked_fly_of_interest;
        
        genotype_mask = strcmp(exp_info.genotype_str, genotype_of_interest{j});
        
        scatter(...
            lunge_interval_all_cell{period_mask, is_attacked_mask, genotype_mask}./fps, ...
            featAll_cell{feat_mask, stat_mask, period_mask, is_attacked_mask, genotype_mask}, ...
            12, marker_opts{j}, 'LineWidth', 1);
        legends = [legends, sprintf('%s in %s-v-%s', genotype, exp_info.genotype_str{1}, exp_info.genotype_str{2})];
        
        ax = gca;
        ax.ColorOrderIndex = j;

        lunge_freq = length(lunge_start_all_cell{period_mask, is_attacked_mask, genotype_mask})/(length_of_exp*60*num_of_flies(genotype_mask));
        fprintf('Lunge frequency %f per sec \n', lunge_freq);
        ylims = ylim; 
        plot((1/lunge_freq).*ones(100, 1), linspace(ylims(1), ylims(2), 100), 'LineStyle', '--');
        legends = [legends, sprintf('theoretical interval of evenly-spaced lunges at %.2f Hz', lunge_freq)];
    end
    
    hold off; 
    set(gca, 'XScale', 'log');
    xlabel('Interval between two lunges (in sec)');
    ylim([0, pi+0.01]);
    if strcmp(stat_of_interest{i}, 'timepoint')
        ylabel(sprintf('%s at %.1f sec after lunge', feature_of_interest{i}, abs(max(period_of_interest))./fps), 'Interpreter', 'none');
        ylims = ylim; 
        line(0.5.*ones(100, 1), linspace(ylims(1), ylims(2), 100), 'LineStyle', '--', 'Color', 'k');
        legends = [legends, sprintf('cutoff on interval between lunge at %.1f sec', abs(max(period_of_interest))./fps)];
    else
        ylabel(strcat(feature_of_interest{i}, '_', stat_of_interest{i}), 'Interpreter', 'none');
        ylims = ylim; 
        line((1/fps).*ones(100, 1), linspace(ylims(1), ylims(2), 100), 'LineStyle', '--', 'Color', 'k');
        legends = [legends, sprintf('shortest resolvable interval at recording rate %dHz', fps)];
    end
    legend(legends, 'Interpreter', 'none');
end
