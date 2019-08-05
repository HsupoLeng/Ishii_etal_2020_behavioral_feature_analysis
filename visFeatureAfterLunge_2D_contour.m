close all; 
load('featAll_and_feat_probs_dist_remove_outliers.mat', ...
    'feature_options', 'stat_options', 'period_options', 'is_attacked_fly', 'featAll_cell', 'period_mask_all_cell');

% Visualize wildtype and mutant feature statistic distributions using heat
% map based on kernel density estimate from the statistic data points
% Compare by visual examination
period_options = fliplr(period_options);
plot_periods = [10, 30, 60];

chosen_is_attacked_fly = true; 
setup_idx_tuples = [repelem((1:length(feature_options))', length(stat_options)), ...
    repmat((1:length(stat_options))', length(feature_options), 1)];

% For now, we only look at combination of different features
setup_combs = nchoosek(1:size(setup_idx_tuples, 1), 2);
comb_mask = diff(floor((setup_combs-1)./length(stat_options)), 1, 2) ~= 0;
setup_combs = setup_combs(comb_mask, :); 

chosen_is_attacked_fly_idx = find(is_attacked_fly == chosen_is_attacked_fly);
ks_num_points = 90; 

% Unify density estimate range for all plots
ks_support_global_min_cell = cell(length(feature_options), length(feature_options));
ks_support_global_max_cell = cell(length(feature_options), length(feature_options));
caxis_common_cell = cell(length(feature_options), length(feature_options));
for i=1:length(feature_options)
    for j=1:length(feature_options)
        ks_support_global_min_cell{i, j} = [inf, inf];
        ks_support_global_max_cell{i, j} = [0, 0];
        caxis_common_cell{i, j} = [inf, 0];
    end
end

for l=1:length(period_options)
    chosen_period = period_options(l);
    if ~ismember(chosen_period, plot_periods)
        continue;
    end

    for i=1:size(setup_combs, 1) % Iterate over different combinations of statistic
        setup_comb = setup_combs(i, :);
        % do not include var and delta in the unification process
        if any(ismember(setup_comb, [5:5:20, 4:5:19]))
            continue; 
        end
        
        feat_pair_mat_cell = cell(2, 1);
        for j=1:2 % Iterate over wildtype and mutant
            feat_pair_mat = [];
            for k=1:length(setup_comb)
                setup_tuple = setup_idx_tuples(setup_comb(k), :);
                feat_pair_mat(:, k) = featAll_cell{setup_tuple(1), setup_tuple(2), length(period_options) - l + 1, chosen_is_attacked_fly_idx, j};    
            end
%             if remove_short_ili_seq
%                 feat_pair_mat = feat_pair_mat(period_mask_all_cell{l, chosen_is_attacked_fly_idx, j}, :);
%             end
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
        
%         [ks_support_mesh_x, ks_support_mesh_y] = meshgrid(linspace(ks_support_global_min(1), ks_support_global_max(1), ks_num_points), ...
%             linspace(ks_support_global_min(2), ks_support_global_max(2), ks_num_points)); 
%         
%         caxis_common = caxis_common_cell{setup_idx_tuples(setup_comb(1), 1), setup_idx_tuples(setup_comb(2), 1)};
%         % Run kernel density estimate
%         for j=1:2
%             [f, xi] = ksdensity(feat_pair_mat_cell{j}, [ks_support_mesh_x(:), ks_support_mesh_y(:)]);
%             if min(f) < caxis_common(1)
%                 caxis_common(1) = min(f);
%             end
%             if max(f) > caxis_common(2)
%                 caxis_common(2) = max(f);
%             end
%         end
%         caxis_common_cell{setup_idx_tuples(setup_comb(1), 1), setup_idx_tuples(setup_comb(2), 1)} = caxis_common;
    end
end

% Draw the plots
c_map = parula(length(plot_periods));
c_map = flipud(c_map);
for i=1:size(setup_combs, 1) % Iterate over different combinations of statistic
    setup_comb = setup_combs(i, :);
    % do not include var and delta in the unification process
    if any(ismember(setup_comb, [5:5:20, 4:5:19]))
        continue; 
    end
    
    ax_objs = cell(2, 1);
    plot_period_accum = 0;
    for l=1:length(period_options)
        chosen_period = period_options(l);
        if ~ismember(chosen_period, plot_periods)
            continue;
        else
            plot_period_accum = plot_period_accum + 1; 
        end
        
        feat_pair_mat_cell = cell(2, 1);
        for j=1:2 % Iterate over wildtype and mutant
            feat_pair_mat = [];
            for k=1:length(setup_comb)
                setup_tuple = setup_idx_tuples(setup_comb(k), :);
                feat_pair_mat(:, k) = featAll_cell{setup_tuple(1), setup_tuple(2), length(period_options) - l + 1, chosen_is_attacked_fly_idx, j};    
            end
            feat_pair_mat_cell{j} = feat_pair_mat;
        end
        
        ks_support_global_min = ks_support_global_min_cell{setup_idx_tuples(setup_comb(1), 1), setup_idx_tuples(setup_comb(2), 1)};
        ks_support_global_max = ks_support_global_max_cell{setup_idx_tuples(setup_comb(1), 1), setup_idx_tuples(setup_comb(2), 1)};
        
        if any(ismember(setup_comb, 1:5))
            ks_support_global_max(1) = 10;
        end
        if any(ismember(setup_comb, 15:20))
            ks_support_global_max(2) = 20;
        end
        [ks_support_mesh_x, ks_support_mesh_y] = meshgrid(linspace(ks_support_global_min(1), ks_support_global_max(1), ks_num_points), ...
            linspace(ks_support_global_min(2), ks_support_global_max(2), ks_num_points)); 

        fig_obj = figure(i);
        for j=1:2  
            [f, xi] = ksdensity(feat_pair_mat_cell{j}, [ks_support_mesh_x(:), ks_support_mesh_y(:)]);
            if isempty(ax_objs{j})
                ax_objs{j} = subplot(1,2,j);
                [~, c_obj] = contour(ax_objs{j}, reshape(xi(:,1), ks_num_points, ks_num_points), reshape(xi(:,2), ks_num_points, ks_num_points), ...
                    reshape(f, ks_num_points, ks_num_points), [0.03, 0.01], 'b', 'LineWidth', 2); drawnow;
                set(c_obj, 'LineColor', uint8(c_map(plot_period_accum, :)'.*255));
            else
                hold(ax_objs{j}, 'on'); 
                if plot_period_accum == length(plot_periods)
                    contour_txt_opt = 'on';
                else
                    contour_txt_opt = 'off';
                end
                [~, c_obj] = contour(ax_objs{j}, reshape(xi(:,1), ks_num_points, ks_num_points), reshape(xi(:,2), ks_num_points, ks_num_points), ...
                    reshape(f, ks_num_points, ks_num_points), [0.03, 0.01], 'b', 'LineWidth', 2, 'ShowText', contour_txt_opt); drawnow; 
                set(c_obj, 'LineColor', uint8(c_map(plot_period_accum, :)'.*255));
                hold(ax_objs{j}, 'off'); 
            end
            
%             surf(reshape(xi(:,1), ks_num_points, ks_num_points), reshape(xi(:,2), ks_num_points, ks_num_points), reshape(f, ks_num_points, ks_num_points), 'EdgeColor', 'none');
%             colormap(teals); 
%             view(2);
             

            axis tight;
            for k = 1:length(setup_comb)
                setup_tuple = setup_idx_tuples(setup_comb(k), :);
                if k == 1
                    xlabel(strrep(strcat(feature_options{setup_tuple(1)}, '-', stat_options{setup_tuple(2)}), '_', '-'));
                    if contains(feature_options{setup_tuple(1)}, 'angle')
                        xticks((0:4).*(pi/4));
                        xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
                        xlim([0, pi]);
                    end
                else
                    ylabel(strrep(strcat(feature_options{setup_tuple(1)}, '-', stat_options{setup_tuple(2)}), '_', '-'));
                    if contains(feature_options{setup_tuple(1)}, 'angle')
                        yticks((0:4).*(pi/4));
                        yticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
                        ylim([0, pi]);
                    end
                end
            end

            if j == 1
                title('Wild-type');
            else
                title('Mutant');
            end
        end

%         caxis_common = caxis_common_cell{setup_idx_tuples(setup_comb(1), 1), setup_idx_tuples(setup_comb(2), 1)};
%         for j=1:2
%             caxis(ax_objs{j}, caxis_common);
%             if j == 2
%                 original_pos = get(ax_objs{j}, 'Position');
%                 colorbar; 
%                 set(ax_objs{j}, 'Position', original_pos);
%             end
%         end
%         % linkaxes([ax_objs{:}], 'xy');
%         set(fig_obj, 'Position', get(0, 'Screensize'));
%         suptitle(sprintf('Period %d', chosen_period));

        if chosen_is_attacked_fly
            fly_str = 'attacked';
        else
            fly_str = 'attacking';
        end
        
    end
    
    for j=1:2
        lgd_obj = legend(ax_objs{j}); 
        temp_p_objs = lgd_obj.PlotChildren; 
        legend('off');
        legend(ax_objs{j}, flipud(temp_p_objs), cellstr(num2str(flipud(plot_periods)')));
    end
    saveas(double(fig_obj), sprintf('contour-wt_vs_mutant-%s_vs_%s-%s_fly.png', ...
            get(get(gca, 'ylabel'), 'String'), get(get(gca, 'xlabel'), 'String'), fly_str));
    saveas(double(fig_obj), sprintf('contour_heat_map-wt_vs_mutant-%s_vs_%s-%s_fly.eps', ...
        get(get(gca, 'ylabel'), 'String'), get(get(gca, 'xlabel'), 'String'), fly_str));
end