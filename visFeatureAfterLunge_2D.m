close all; 
load('flymat_info_struct.mat', 'flymat_info_struct');
plot_periods = [-30, -1, 30];

feat_pair_mat_cell_all = cell(length(flymat_info_struct), length(plot_periods), 3, 2);
skip_plot = true; 
for s=1:length(flymat_info_struct)
    flymat_info = flymat_info_struct(s);
    is_same_genotype_pairing = flymat_info.is_same_genotype_pairing;
    % featAll_mat_names will always be of length 1 currently
    if is_same_genotype_pairing
        same_genotype_str = 'same_genotype';
        featAll_mat_names = strcat('featAll_and_feat_probs_dist_remove_outliers-', same_genotype_str, '-', strjoin(flymat_info.flymat, '_'));
        featAll_mat_names = {featAll_mat_names};
    else
        same_genotype_str = 'diff_genotype';
        featAll_mat_names = cellfun(@(flymat_name) strcat('featAll_and_feat_probs_dist_remove_outliers-', same_genotype_str, '-', flymat_name), flymat_info.flymat, 'UniformOutput', false);
    end
    
    
%     featAll_mat_names = {'featAll_and_feat_probs_dist_remove_outliers-FLYMAT_CSMH_SM_vsWT_GM_HeisenbergChamber_May2019', ...
%         'featAll_and_feat_probs_dist_remove_outliers-FLYMAT_CG3385_KO_CS-bc-11_GM_vsWT_GM_HeisenbergChamber_May2019'};
    % featAll_mat_names = {'featAll_and_feat_probs_dist_remove_outliers-FLYMAT_CG3385_KO_CS-bc-11_GM_vsWT_GM_HeisenbergChamber_May2019'};

    load(featAll_mat_names{1}, ...
            'feature_options', 'stat_options', 'period_options', 'is_attacked_fly');
    % Visualize wildtype and mutant feature statistic distributions using heat
    % map based on kernel density estimate from the statistic data points
    % Compare by visual examination
    
    % plot_periods(plot_periods == 0) = [];
    chosen_is_attacked_fly = false; 
    chosen_is_attacked_fly_idx = find(is_attacked_fly == chosen_is_attacked_fly);
    ks_num_points = 90; 
    setup_idx_tuples = [repelem((1:length(feature_options))', length(stat_options)), ...
        repmat((1:length(stat_options))', length(feature_options), 1)];
    % For now, we only look at combination of different features
    % setup_combs = nchoosek(1:size(setup_idx_tuples, 1), 2);
    % comb_mask = diff(floor((setup_combs-1)./length(stat_options)), 1, 2) ~= 0;
    % setup_combs = setup_combs(comb_mask, :); 
    if chosen_is_attacked_fly
        fly_str = 'attacked';
    else
        fly_str = 'attacking';
    end

    featAll_all_experiments_struct = struct();
    for p=1:length(featAll_mat_names)
        featAll_temp_struct = load(featAll_mat_names{p}, 'featAll_cell', 'genotype_str', 'lunge_start_all_cell');
        featAll_all_experiments_struct(p).featAll_cell = featAll_temp_struct.featAll_cell;
        featAll_all_experiments_struct(p).genotype_str = featAll_temp_struct.genotype_str; 
        featAll_all_experiments_struct(p).lunge_start_all_cell = featAll_temp_struct.lunge_start_all_cell; 
    end


    ks_support_global_min_cell = cell(length(feature_options), length(feature_options));
    ks_support_global_max_cell = cell(length(feature_options), length(feature_options));
    caxis_common_cell = cell(length(featAll_all_experiments_struct), length(feature_options), length(feature_options), 2);
    for i=1:length(feature_options)
        for j=1:length(feature_options)
            ks_support_global_min_cell{i, j} = [inf, inf];
            ks_support_global_max_cell{i, j} = [0, 0];
        end
    end

    for p=1:length(featAll_all_experiments_struct)
        for i=1:length(feature_options)
            for j=1:length(feature_options)
                for r=1:2
                    caxis_common_cell{p, i, j, r} = [inf, 0];
                end
            end
        end
    end

    % Unify density estimate range for all plots, across all experiments
    for p=1:length(featAll_mat_names)
        featAll_cell = featAll_all_experiments_struct(p).featAll_cell; 
        genotype_str = featAll_all_experiments_struct(p).genotype_str; 

        for l=1:length(period_options)
            chosen_period = period_options(l);
            if ~ismember(chosen_period, plot_periods)
                continue; 
            end
            if chosen_period < 0 
                setup_combs = [22, 27; 27, 2; 22, 2];
            %     setup_combs = [2, 27];
            %     setup_combs = [31, 36];
            else
                setup_combs = [23, 28; 28, 3; 23, 3];
            %     setup_combs = [3, 28];
            %     setup_combs = [31, 36];
            end

            for i=1:size(setup_combs, 1) % Iterate over different combinations of statistic
                setup_comb = setup_combs(i, :);
                % do not include var and delta in the unification process
                if any(ismember(setup_comb, [5:5:20, 4:5:19]))
                    continue; 
                end

                feat_pair_mat_cell = cell(2, 1);
                for j=1:2 % Iterate over wildtype and mutant
                    if strcmp(genotype_str{j}, 'WT_GM')
                        continue; 
                    end

                    feat_pair_mat = [];
                    for k=1:length(setup_comb)
                        setup_tuple = setup_idx_tuples(setup_comb(k), :);
                        feat_pair_mat(:, k) = featAll_cell{setup_tuple(1), setup_tuple(2), l, chosen_is_attacked_fly_idx, j};    
                    end
                    feat_pair_mat_cell{j} = feat_pair_mat;
                end

                ks_support_min = min(vertcat(feat_pair_mat_cell{:}));
                ks_support_min = (1-0.2.*(sign(ks_support_min))).*ks_support_min;
                ks_support_max = max(vertcat(feat_pair_mat_cell{:}));
                ks_support_max = (1+0.2.*(sign(ks_support_max))).*ks_support_max;

                ks_support_global_min = ks_support_global_min_cell{setup_idx_tuples(setup_comb(1), 1), setup_idx_tuples(setup_comb(2), 1)};
                ks_support_global_max = ks_support_global_max_cell{setup_idx_tuples(setup_comb(1), 1), setup_idx_tuples(setup_comb(2), 1)};

                if any(ks_support_min < ks_support_global_min)
                    ks_support_global_min(ks_support_min < ks_support_global_min) = ...
                        ks_support_min(ks_support_min < ks_support_global_min);
                    ks_support_global_min_cell{setup_idx_tuples(setup_comb(1), 1), setup_idx_tuples(setup_comb(2), 1)} = ...
                        ks_support_global_min; 
                end
                if any(ks_support_max > ks_support_global_max)
                    ks_support_global_max(ks_support_max > ks_support_global_max) = ...
                        ks_support_max(ks_support_max > ks_support_global_max);
                    ks_support_global_max_cell{setup_idx_tuples(setup_comb(1), 1), setup_idx_tuples(setup_comb(2), 1)} = ...
                        ks_support_global_max; 
                end

                if any(ismember(setup_comb, 1:5))
                    ks_support_global_min(2) = 1;
                    ks_support_global_max(2) = 10;
                end
                if all(ismember(setup_comb, 21:30))
                    ks_support_global_min(1) = 0; 
                    ks_support_global_min(2) = 0;
                end
                if any(ismember(setup_comb, 15:20)) || any(ismember(setup_comb, 31:35))
                    ks_support_global_max(1) = 20;
                    ks_support_global_max(2) = 20;
                end
                [ks_support_mesh_x, ks_support_mesh_y] = meshgrid(linspace(ks_support_global_min(1), ks_support_global_max(1), ks_num_points), ...
                    linspace(ks_support_global_min(2), ks_support_global_max(2), ks_num_points)); 

                caxis_common = caxis_common_cell{p, setup_idx_tuples(setup_comb(1), 1), setup_idx_tuples(setup_comb(2), 1), (chosen_period < 0) + 1};
                % Run kernel density estimate
                for j=1:2
                    if strcmp(genotype_str{j}, 'WT_GM')
                        continue; 
                    end

                    [f, xi] = ksdensity(feat_pair_mat_cell{j}, [ks_support_mesh_x(:), ks_support_mesh_y(:)]);
                    if min(f) < caxis_common(1)
                        caxis_common(1) = min(f);
                    end
                    if max(f) > caxis_common(2)
                        caxis_common(2) = max(f);
                    end
                end
                caxis_common_cell{p, setup_idx_tuples(setup_comb(1), 1), setup_idx_tuples(setup_comb(2), 1), (chosen_period < 0) + 1} = caxis_common;
            end
        end
    end

    % Draw the plots
    % c_map = parula(length(plot_periods));
    for p=1:length(featAll_mat_names)
        featAll_cell = featAll_all_experiments_struct(p).featAll_cell; 
        genotype_str = featAll_all_experiments_struct(p).genotype_str; 
        lunge_start_all_cell = featAll_all_experiments_struct(p).lunge_start_all_cell; 
        quadrant_count_cell = cell(length(period_options), size(setup_combs, 1), 2);

        for l=1:length(period_options)
            chosen_period = period_options(l);
            if ~ismember(chosen_period, plot_periods)
                continue;
            end
            if chosen_period < 0 
                setup_combs = [22, 27; 27, 2; 22, 2];
            %     setup_combs = [2, 27];
            %     setup_combs = [31, 36];
            else
                setup_combs = [23, 28; 28, 3; 23, 3];
            %     setup_combs = [3, 28];
            %     setup_combs = [31, 36];
            end
                 
            for i=1:size(setup_combs, 1) % Iterate over different combinations of statistic
                setup_comb = setup_combs(i, :);
                % do not include var and delta in the unification process
                if any(ismember(setup_comb, [5:5:20, 4:5:19]))
                    continue; 
                end

                ax_objs = cell(2, 1);

                ks_support_global_min = ks_support_global_min_cell{setup_idx_tuples(setup_comb(1), 1), setup_idx_tuples(setup_comb(2), 1)};
                ks_support_global_max = ks_support_global_max_cell{setup_idx_tuples(setup_comb(1), 1), setup_idx_tuples(setup_comb(2), 1)};

                if any(ismember(setup_comb, 1:5))                
                    ks_support_global_min(2) = 1;
                    ks_support_global_max(2) = 10;
                end
                if any(ismember(setup_comb, 21:30))
                    ks_support_global_min(ismember(setup_comb, 21:30)) = 0; 
                end
                if any(ismember(setup_comb, 15:20)) || any(ismember(setup_comb, 31:35))
                    ks_support_global_min(1) = 0;
                    ks_support_global_min(2) = 0;
                    ks_support_global_max(1) = 20;
                    ks_support_global_max(2) = 20;
                end
                [ks_support_mesh_x, ks_support_mesh_y] = meshgrid(linspace(ks_support_global_min(1), ks_support_global_max(1), ks_num_points), ...
                    linspace(ks_support_global_min(2), ks_support_global_max(2), ks_num_points)); 

                if ismember(2, setup_comb) || ismember(3, setup_comb)
                    xaxis_opt = 'angle';
                    yaxis_opt = 'dist'; 
                    yaxis_scale_opt = 10 - ks_support_global_min(2);
                    xaxis_scale_opt = pi;
                    xaxis_draw_offset = ks_support_global_min(1);
                    yaxis_draw_offset = ks_support_global_min(2);
                elseif ismember(31, setup_comb)
                    xaxis_opt = 'vel';
                    yaxis_opt = 'vel';
                    xaxis_scale_opt = 20;
                    yaxis_scale_opt = 20;
                    xaxis_draw_offset = 0;
                    yaxis_draw_offset = 0;
                else
                    xaxis_opt = 'angle'; 
                    yaxis_opt = 'angle';
                    xaxis_scale_opt = pi;
                    yaxis_scale_opt = pi;
                    xaxis_draw_offset = ks_support_global_min(1);
                    yaxis_draw_offset = ks_support_global_min(2);
                end
                role_opt = cell(size(setup_comb));
                if any(ismember(setup_comb, [22, 23]))
                    role_opt{ismember(setup_comb, [22, 23])} = 'attacking';
                end
                if any(ismember(setup_comb, [27, 28]))
                    role_opt{ismember(setup_comb, [27, 28])} = 'attacked';
                end

                feat_pair_mat_cell = cell(2, 1);
                for j=1:2 % Iterate over wildtype and mutant
                    if strcmp(genotype_str{j}, 'WT_GM')
                        continue; 
                    end

                    feat_pair_mat = [];
                    for k=1:length(setup_comb)
                        setup_tuple = setup_idx_tuples(setup_comb(k), :);
                        feat_pair_mat(:, k) = featAll_cell{setup_tuple(1), setup_tuple(2), l, chosen_is_attacked_fly_idx, j};    
                    end
                    quadrant_count = zeros(3, 3);
                    quadrant_members_start_frame = cell(3, 3, 2);
                    for m=1:size(quadrant_count, 1)
                        for n=1:size(quadrant_count, 2)
                            quadrant_mask = bitand(bitand(feat_pair_mat(:, 1)< xaxis_draw_offset + m*xaxis_scale_opt/size(quadrant_count, 1), feat_pair_mat(:, 1)> xaxis_draw_offset + (m-1)*xaxis_scale_opt/size(quadrant_count, 1)), ...
                                bitand(feat_pair_mat(:, 2) < yaxis_draw_offset + n*yaxis_scale_opt/size(quadrant_count, 2), feat_pair_mat(:, 2) > yaxis_draw_offset + (n-1)*yaxis_scale_opt/size(quadrant_count, 2)));
                            quadrant_members_start_frame{m, n, j} = lunge_start_all_cell{l, 1, j}(quadrant_mask);
                            quadrant_count(m, n) = sum(quadrant_mask);                   
                        end
                    end
                    quadrant_count_cell{l, i, j} = quadrant_count;
                    feat_pair_mat_cell{j} = feat_pair_mat; 
                    feat_pair_mat_cell_all{s, l, i, j} = feat_pair_mat;
                end
                
                if skip_plot
                    continue; 
                end
                
                fig_obj = figure(i);
                for j=1:2  
                    if strcmp(genotype_str{j}, 'WT_GM')
                        continue; 
                    end

                    [f, xi] = ksdensity(feat_pair_mat_cell{j}, [ks_support_mesh_x(:), ks_support_mesh_y(:)]);

    %                 ax_objs{j} = subplot(1, 2, j);
                    ax_objs{j}  = gca;
                    surf(reshape(xi(:,1), ks_num_points, ks_num_points), reshape(xi(:,2), ks_num_points, ks_num_points), reshape(f, ks_num_points, ks_num_points), 'EdgeColor', 'none');
                    colormap(teals); 
                    view(2);
                    axis square;
                    axis tight;

                    for k = 1:length(setup_comb)
                        setup_tuple = setup_idx_tuples(setup_comb(k), :);
                        if k == 1
                            if strcmp(xaxis_opt, 'angle')
                                xlabel(sprintf('Facing angle of %s fly', role_opt{k}));
                            elseif strcmp(xaxis_opt, 'vel')
                                xlabel(sprintf('Mean velocity of %s fly', role_opt{k}));
                            elseif strcmp(xaxis_opt, 'dist')
                                xlabel('Distance between two flies');
                            end
                            if contains(feature_options{setup_tuple(1)}, 'angle')
                                xticks((0:3).*(pi/3));
                                xticklabels({'0', '\pi/3', '2\pi/3', '\pi'});
                                xlim([0, pi]);
                            else
                                xticks(linspace(ks_support_global_min(1), ks_support_global_max(1), 4));
                                xtickformat('%.1f');
                            end
                        else
                            if strcmp(yaxis_opt, 'angle')
                                ylabel(sprintf('Facing angle of %s fly', role_opt{k}));
                            elseif strcmp(yaxis_opt, 'vel')
                                ylabel(sprintf('Mean velocity of %s fly', role_opt{k}));
                            elseif strcmp(yaxis_opt, 'dist')
                                ylabel('Distance between two flies');
                            end

                            if contains(feature_options{setup_tuple(1)}, 'angle')
                                yticks((0:3).*(pi/3));
                                yticklabels({'0', '\pi/3', '2\pi/3', '\pi'});
                                ylim([0, pi]);
                            else
                                yticks(linspace(ks_support_global_min(2), ks_support_global_max(2), 4));
                                ytickformat('%.1f');
                            end
                        end
                    end

                    quadrant_count = quadrant_count_cell{l, i, j};
                    for m = 1:size(quadrant_count, 1)
                        for n = 1:size(quadrant_count, 2)
                            text(xaxis_draw_offset + m*xaxis_scale_opt/size(quadrant_count, 1) - 0.5*xaxis_scale_opt/size(quadrant_count, 1), yaxis_draw_offset + n*yaxis_scale_opt/size(quadrant_count, 2) - 0.5*yaxis_scale_opt/size(quadrant_count, 1), 1, ...
                                {sprintf('  %d%%', round(100*quadrant_count(m, n)/sum(sum(quadrant_count)))), sprintf('N = %d', quadrant_count(m, n))});
                        end
                        if m ~= size(quadrant_count, 1)
                            line([xaxis_draw_offset, xaxis_draw_offset + xaxis_scale_opt], [yaxis_draw_offset + m*yaxis_scale_opt/size(quadrant_count, 1), yaxis_draw_offset + m*yaxis_scale_opt/size(quadrant_count, 1)], [1, 1], 'LineStyle', '--');
                            line([xaxis_draw_offset + m*xaxis_scale_opt/size(quadrant_count, 1), xaxis_draw_offset + m*xaxis_scale_opt/size(quadrant_count, 1)], [yaxis_draw_offset, yaxis_draw_offset + yaxis_scale_opt], [1, 1], 'LineStyle', '--');
                        end
                    end

                    title(strrep(genotype_str{j}, '_', '-'));

                    if is_same_genotype_pairing
                        pairing_str = strjoin(repmat(genotype_str(j), 2, 1), '_vs_');
                    else
                        pairing_str = strjoin(genotype_str, '_vs_');
                    end

%                     original_pos = get(ax_objs{j}, 'Position');
%                     set(ax_objs{j}, 'Position', original_pos);
                    caxis_common = caxis_common_cell{p, setup_idx_tuples(setup_comb(1), 1), setup_idx_tuples(setup_comb(2), 1), (chosen_period<0) + 1};
                    caxis(ax_objs{j}, caxis_common);
                    colorbar;

                    set(fig_obj, 'Position', get(0, 'Screensize'));
                    if chosen_period < 0
                        period_sign = 'm';
                    else
                        period_sign = 'p';
                    end
                    set(gcf,'renderer','Painters');
                    saveas(double(fig_obj), sprintf('teal_heat_map-%s-%s-%s_vs_%s-period_%s%d-%s_fly.png', ...
                        pairing_str, genotype_str{j}, get(get(gca, 'ylabel'), 'String'), get(get(gca, 'xlabel'), 'String'), period_sign, abs(chosen_period), fly_str));
                    saveas(double(fig_obj), sprintf('teal_heat_map-%s-%s-%s_vs_%s-period_%s%d-%s_fly.eps', ...
                        pairing_str, genotype_str{j}, get(get(gca, 'ylabel'), 'String'), get(get(gca, 'xlabel'), 'String'), period_sign, abs(chosen_period), fly_str), 'epsc');

                    quadrant_members_start_frame_tester = quadrant_members_start_frame(:, :, 1);
                    save(sprintf('quadrant_members_start_frame-%s-%s-%s_vs_%s-period_%d-%s_fly.mat', ...
                            pairing_str, genotype_str{j}, get(get(gca, 'ylabel'), 'String'), get(get(gca, 'xlabel'), 'String'), chosen_period, fly_str), 'quadrant_members_start_frame_tester');
                end
            end         
        end
    end
end

%% Compute 1D and 2D distribution p-values
pairs_to_test = {{[1, 1], [1, 2]}, {[2, 1], [2, 2]}, {[3, 2], [4, 2]}};
feature_pair_str = {'facing angle of self vs. facing angle of other', 'facing angle of other vs. distance', 'facing angle of self vs. distance'};
kstest2d_p_val_all = zeros(length(pairs_to_test), length(plot_periods), length(feature_pair_str));
ranksum_p_val_all = zeros(length(pairs_to_test), length(plot_periods), length(feature_pair_str), 2);
kstest1d_p_val_all = zeros(size(ranksum_p_val_all));
for s=1:length(pairs_to_test)
    pair = pairs_to_test{s};
    for l=1:size(kstest2d_p_val_all, 2)
        for i=1:size(kstest2d_p_val_all, 3)
            feat_pair_mat_1 = feat_pair_mat_cell_all{pair{1}(1), l, i, pair{1}(2)}; 
            feat_pair_mat_2 = feat_pair_mat_cell_all{pair{2}(1), l, i, pair{2}(2)}; 
            smpl_sz = min(size(feat_pair_mat_1, 1), size(feat_pair_mat_2, 1));
%             [kstest2d_p_val_all(s, l, i), ~] = ...  
%                 kstest2d(feat_pair_mat_1(randi(size(feat_pair_mat_1, 1), smpl_sz, 1), :), ...
%                 feat_pair_mat_2(randi(size(feat_pair_mat_2, 1), smpl_sz, 1), :));
            
            % Use large (> 300) but unequal number of points in each sample 
            [kstest2d_p_val_all(s, l, i), ~] = ...  
                kstest2d(feat_pair_mat_1, feat_pair_mat_2); % 2D Kolmogorov-Smirnov
            for j=1:2
                ranksum_p_val_all(s, l, i, j) = ranksum(feat_pair_mat_1(:, j), feat_pair_mat_2(:, j));
            end
            for j=1:2
                [~, kstest1d_p_val_all(s, l, i, j)] = kstest2(feat_pair_mat_1(:, j), feat_pair_mat_2(:, j));
            end
            pairing_type_str_cell = cell(size(pair));
            for k=1:length(pairing_type_str_cell)
                if flymat_info_struct(pair{k}(1)).is_same_genotype_pairing
                    pairing_type_str_cell{k} = sprintf('(%s vs. %s)', ...
                        flymat_info_struct(pair{k}(1)).genotype_str{pair{k}(2)}, ...
                        flymat_info_struct(pair{k}(1)).genotype_str{pair{k}(2)});
                else
                    pairing_type_str_cell{k} = sprintf('(%s vs. %s)', ...
                        flymat_info_struct(pair{k}(1)).genotype_str{1}, ...
                        flymat_info_struct(pair{k}(1)).genotype_str{2});
                end
            end
            fprintf('===== =====\n');
            fprintf('%s \n Genotype %s in %s and %s in %s\n Frame(s) w.r.t. lunge %d: p-value %.3e \n', ...
                feature_pair_str{i}, ...
                flymat_info_struct(pair{1}(1)).genotype_str{pair{1}(2)}, ...
                pairing_type_str_cell{1},...
                flymat_info_struct(pair{2}(1)).genotype_str{pair{2}(2)}, ...
                pairing_type_str_cell{2}, ...
                plot_periods(l), kstest2d_p_val_all(s, l, i));
            feature_pair = regexp(feature_pair_str{i}, 'vs. ', 'split');
            for j=1:2
                fprintf('%s 1D test: ranksum p-val %.3e; KS p-val %.3e\n', ...
                    feature_pair{j}, ...
                    ranksum_p_val_all(s, l, i, j), kstest1d_p_val_all(s, l, i, j));
            end
        end
    end
end
save('stat_test_p_vals.mat', 'ranksum_p_val_all', 'kstest1d_p_val_all', 'kstest2d_p_val_all', 'pairs_to_test', 'feature_pair_str');
