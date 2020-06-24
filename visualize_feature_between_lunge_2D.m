% Visualize facing angle vs. lunge interval 
% Taking data compiled using extract_features_and_distance, stored in
% flymat_info_struct
load('flymat_info_struct');
exp_of_interest = [3, 4];
genotype_of_interest = {'WT_SM', 'CG3385_KO_CS-bc-11_GM'};
feature_of_interest = {'facing_angle_in_lunge_interval'};
stat_of_interest = {'max'};
period_of_interest = [1, 30];
is_attacked_fly_of_interest = false; 
marker_opts = {'^', '^'};
fps = 60; 
length_of_exp = 30; 
num_of_bins = 3; 
interval_scales = 10.^(-2:3)';
interval_scales = [interval_scales(1:end-1), interval_scales(2:end)];
interval_scale_for_bin = [2, 3, 3];
uniform_scale = 10; 
large_scale = true; 
interval_scale_for_feature = [0, 2, 20, 1000];
feature_num_of_bins = 12; 
zoom_in = true; 
skip_type = 1; 

interval_binned_by_feature_all = cell(num_of_bins, length(exp_of_interest), length(feature_of_interest));
feature_binned_by_interval_all = cell(num_of_bins, length(exp_of_interest), length(feature_of_interest));

for i=1:length(feature_of_interest)
    s = figure(); 
    h = figure();
    hf = figure(); 
    p = figure();
    
    s_legends = {}; 
    h_legends = {};
    for j=1:length(exp_of_interest)
        exp_info = flymat_info_struct(exp_of_interest(j));
        genotype = genotype_of_interest{j};
        
        if exp_info.is_same_genotype_pairing
            same_genotype_str = 'same_genotype';
        else
            same_genotype_str = 'diff_genotype';
        end
        load(strcat('featAll_and_feat_probs_dist_keep_outliers-', ...
            same_genotype_str, '-', strjoin(exp_info.flymat, '_'), '.mat'), ...
            'featAll_cell', 'lunge_interval_all_cell', ...
            'feature_options', 'stat_options', 'period_options', 'is_attacked_fly', ...
            'lunge_start_all_cell', 'num_of_flies');

        feat_mask = strcmp(feature_options, feature_of_interest{i});
        stat_mask = strcmp(stat_options, stat_of_interest{i});
        period_mask = cellfun(@(period_pair) isequal(period_pair, period_of_interest), period_options); 
        is_attacked_mask = is_attacked_fly == is_attacked_fly_of_interest;
        genotype_mask = strcmp(exp_info.genotype_str, genotype_of_interest{j});
        
        lunge_interval_all = lunge_interval_all_cell{period_mask, is_attacked_mask, genotype_mask}./fps;
        feat_all = featAll_cell{feat_mask, stat_mask, period_mask, is_attacked_mask, genotype_mask};
        
        % Draw feature vs. lunge interval scatter plot
        figure(s);
        hold on;
        
        ax = gca;
        ax.ColorOrderIndex = j;
        if ~skip_type || ( skip_type && j ~= skip_type) 
            scatter(lunge_interval_all, feat_all, 12, marker_opts{j}, 'LineWidth', 1);
            s_legends = [s_legends, sprintf('%s in %s-v-%s', genotype, exp_info.genotype_str{1}, exp_info.genotype_str{2})];
        end
        
        ax = gca;
        ax.ColorOrderIndex = j;

        lunge_freq = length(lunge_start_all_cell{period_mask, is_attacked_mask, genotype_mask})/(length_of_exp*60*num_of_flies(genotype_mask));
        fprintf('Lunge frequency %f per sec \n', lunge_freq);
        ylims = ylim; 
        if (skip_type && skip_type ~= j) || (~skip_type && j == 2)
            plot((1/lunge_freq).*ones(100, 1), linspace(ylims(1), ylims(2), 100), 'LineStyle', '--');
            s_legends = [s_legends, sprintf('theoretical interval of evenly-spaced lunges at %.2f Hz', lunge_freq)];
        end
        for k=1:num_of_bins
            feat_bin_ratio = histcounts(feat_all, (0:num_of_bins).*pi/num_of_bins, 'Normalization', 'Probability');
        end
        
        set(gca, 'XScale', 'log');
%         axis square; 
        xlabel('Interval between two lunges (in sec)');
        xlim([0.01, 1000]);
        ylim([0, pi+0.01]);
        if strcmp(stat_of_interest{i}, 'timepoint')
            ylabel(sprintf('%s at %.1f sec after lunge', feature_of_interest{i}, abs(max(period_of_interest))./fps), 'Interpreter', 'none');
            ylims = ylim; 
            line(0.5.*ones(100, 1), linspace(ylims(1), ylims(2), 100), 'LineStyle', '--', 'Color', 'k');
            s_legends = [s_legends, sprintf('cutoff on interval between lunge at %.1f sec', abs(max(period_of_interest))./fps)];
        else
            ylabel(strcat(feature_of_interest{i}, '_', stat_of_interest{i}), 'Interpreter', 'none');
            ylims = ylim;
            if (skip_type && skip_type ~= j) || (~skip_type && j == 2)
                line((1/fps).*ones(100, 1), linspace(ylims(1), ylims(2), 100), 'LineStyle', '--', 'Color', 'k');
                s_legends = [s_legends, sprintf('shortest resolvable interval at recording rate %dHz', fps)];
            end
        end
        legend(s_legends, 'Interpreter', 'none');
        if (skip_type && skip_type ~= j) || (~skip_type && j == 2)
            set(gcf, 'Renderer', 'painters');
            if ~skip_type
                skip_str = 'overlaid';
            else
                skip_str = genotype_of_interest{j};
            end
            saveas(gcf, sprintf('facing_angle_max_vs_lunge_interval-%s_triangle.eps', skip_str), 'epsc');
        end
        
        % Draw histogram of interval for samples binned by the feature
        for k=1:num_of_bins
            feat_mask = bitand(feat_all >= (k-1)*pi/num_of_bins, feat_all < k*pi/num_of_bins);
            lunge_interval_selected = lunge_interval_all(feat_mask);
            
            if strcmp(stat_of_interest, 'max')
                if zoom_in
                    edges = linspace(interval_scales(interval_scale_for_bin(k), 1), 2*interval_scales(interval_scale_for_bin(k), 2), 20); 
                    edges = [0, edges]; 
                else
                    edges = [0:uniform_scale:200, 1000];
                end
            elseif strcmp(stat_of_interest, 'timepoint')
                edges = 0:1:20; 
            end
                
            lunge_binned_count = histcounts(lunge_interval_selected(lunge_interval_selected < max(edges)), edges);
            lunge_binned_count = lunge_binned_count ./ length(lunge_interval_all);
            interval_binned_by_feature_all{k, j, i} = lunge_interval_selected; 

            figure(h);
            subplot(num_of_bins, 1, num_of_bins - k + 1);
            if (skip_type && skip_type ~= j) || (~skip_type && j == 2)
                ax = gca;
                ax.ColorOrderIndex = j;
                hold on; 
                b = bar(lunge_binned_count);
                b.FaceAlpha = 0.5;
                h_legends = legend('show');
                h_legends.String = [h_legends.String(1:end-1), sprintf('%s (N = %d / %d)', genotype, length(lunge_interval_selected(lunge_interval_selected < max(edges))), length(lunge_interval_all))];
                h_legends.Interpreter = 'none';
            end
           
            if j == 2
                [~, p_val_unequal] = kstest2(interval_binned_by_feature_all{k, 2, i}, interval_binned_by_feature_all{k, 1, i}, 'Tail', 'unequal');
                [~, p_val_mu_small] = kstest2(interval_binned_by_feature_all{k, 2, i}, interval_binned_by_feature_all{k, 1, i}, 'Tail', 'larger');
                [~, p_val_mu_large] = kstest2(interval_binned_by_feature_all{k, 2, i}, interval_binned_by_feature_all{k, 1, i}, 'Tail', 'smaller');
                title(sprintf('Lunge interval distributions (%s \\in [%d\\pi/%d, %d\\pi/%d]) p-val unequal %.4f, mutant interval smaller %.4f, larger %.4f', ...
                    strcat('facing-angle', '-', stat_of_interest{i}), k-1, num_of_bins, k, num_of_bins, p_val_unequal, p_val_mu_small, p_val_mu_large));
            end

            xticks((1:length(edges))-0.5);
            xticklabels(cellstr(num2str(edges')));
            xlabel('Lunge interval (in sec)');
            ylabel('Probability');
            ylim([0, 0.35]);
            
            if (skip_type && skip_type ~= j) || (~skip_type && j == 2)
                set(gcf, 'Renderer', 'painters');
                if ~skip_type
                    skip_str = 'overlaid';
                else
                    skip_str = genotype_of_interest{j};
                end
            end
        end
        
        % Draw lunge interval-binned feature histogram
        for k=1:num_of_bins
            lunge_interval_mask = bitand(lunge_interval_all >= interval_scale_for_feature(k), ...
                lunge_interval_all < interval_scale_for_feature(k+1));
            feature_selected = feat_all(lunge_interval_mask);
            [feature_binned_count, edges] = histcounts(feature_selected, 12, 'BinLimits', [0, pi]);
            feature_binned_count = feature_binned_count ./ length(feat_all);
            feature_binned_by_interval_all{k, j, i} = feature_selected; 
            
            figure(hf);
            subplot(1, num_of_bins, k);
            if (skip_type && skip_type ~= j) || (~skip_type && j == 2)    
                ax = gca;
                ax.ColorOrderIndex = j;
                hold on; 
                b = bar(feature_binned_count);
                b.FaceAlpha = 0.5;
                h_legends = legend('show');
                h_legends.String = [h_legends.String(1:end-1), sprintf('%s (N = %d / %d)', genotype, length(feature_selected), length(feat_all))];
                h_legends.Interpreter = 'none';
            end
            
            if j == 2
                [~, p_val_unequal] = kstest2(feature_binned_by_interval_all{k, 2, i}, feature_binned_by_interval_all{k, 1, i}, 'Tail', 'unequal');
                [~, p_val_mu_small] = kstest2(feature_binned_by_interval_all{k, 2, i}, feature_binned_by_interval_all{k, 1, i}, 'Tail', 'larger');
                [~, p_val_mu_large] = kstest2(feature_binned_by_interval_all{k, 2, i}, feature_binned_by_interval_all{k, 1, i}, 'Tail', 'smaller');
                title({sprintf('Max facing angle distributions (lunge interval \\in [%d, %d])', ...
                    interval_scale_for_feature(k), interval_scale_for_feature(k+1)), ...
                    sprintf('p-val unequal %.4e', p_val_unequal)});
            end
            
            xticks((1:(feature_num_of_bins+1))-0.5);
            xticklabels(arrayfun(@(n) sprintf('%.2f', n), edges, 'UniformOutput', false));
            xlabel('Facing angle (in radian)');
            xtickangle(-45);
            ylabel('Probability');
            ylim([0, 1])
            
            if (skip_type && skip_type~=j) || (~skip_type && j == 2)
                set(gcf, 'Renderer', 'painters');
            end
        end
    end
end
