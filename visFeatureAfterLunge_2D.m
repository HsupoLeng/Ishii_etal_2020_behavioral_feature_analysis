close all; 
featAll_mat_name = 'featAll_and_feat_probs_dist_remove_outliers-same_genotype_pairing-FLYMAT_CSMH_SM_HeisenbergChamber_Apr-May2019_FLYMAT_CG3385_KO_CS-bc-11_HeisenbergChamber_Apr2019';

load(featAll_mat_name, ...
    'feature_options', 'stat_options', 'period_options', 'is_attacked_fly', 'featAll_cell', 'period_mask_all_cell', 'genotype_str');

% Visualize wildtype and mutant feature statistic distributions using heat
% map based on kernel density estimate from the statistic data points
% Compare by visual examination
plot_periods = [-1, -30];
% plot_periods = [30];
chosen_is_attacked_fly = false; 
setup_idx_tuples = [repelem((1:length(feature_options))', length(stat_options)), ...
    repmat((1:length(stat_options))', length(feature_options), 1)];

% For now, we only look at combination of different features
% setup_combs = nchoosek(1:size(setup_idx_tuples, 1), 2);
% comb_mask = diff(floor((setup_combs-1)./length(stat_options)), 1, 2) ~= 0;
% setup_combs = setup_combs(comb_mask, :); 
if ~any(plot_periods > 0) 
    setup_combs = [22, 27; 2, 27; 2, 22];
%     setup_combs = [2, 27];
elseif ~any(plot_periods < 0)
    setup_combs = [23, 28; 3, 28; 3, 23];
%     setup_combs = [3, 28];
else
    fprintf('Please plot period before lunge and after lunge separately\n');
    return;
end

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

quadrant_count_cell = cell(length(period_options), size(setup_combs, 1), 2);
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
        
        if ismember(2, setup_comb) || ismember(3, setup_comb)
            xaxis_scale_opt = 10; 
        else
            xaxis_scale_opt = pi; 
        end
        
        feat_pair_mat_cell = cell(2, 1);
        for j=1:2 % Iterate over wildtype and mutant
            feat_pair_mat = [];
            for k=1:length(setup_comb)
                setup_tuple = setup_idx_tuples(setup_comb(k), :);
                feat_pair_mat(:, k) = featAll_cell{setup_tuple(1), setup_tuple(2), l, chosen_is_attacked_fly_idx, j};    
            end
            quadrant_count = zeros(3, 3);
            for m=1:size(quadrant_count, 1)
                for n=1:size(quadrant_count, 2)
                    quadrant_count(m, n) = sum(bitand(bitand(feat_pair_mat(:, 1)< m*xaxis_scale_opt/size(quadrant_count, 1), feat_pair_mat(:, 1)> (m-1)*xaxis_scale_opt/size(quadrant_count, 1)), ...
                        bitand(feat_pair_mat(:, 2) < n*pi/size(quadrant_count, 2), feat_pair_mat(:, 2) > (n-1)*pi/size(quadrant_count, 2))));                   
                end
            end
            quadrant_count_cell{l, i, j} = quadrant_count;
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
            ks_support_global_max(1) = 10;
        end
        if any(ismember(setup_comb, 15:20))
            ks_support_global_max(2) = 20;
        end
        [ks_support_mesh_x, ks_support_mesh_y] = meshgrid(linspace(ks_support_global_min(1), ks_support_global_max(1), ks_num_points), ...
            linspace(ks_support_global_min(2), ks_support_global_max(2), ks_num_points)); 
        
        caxis_common = caxis_common_cell{setup_idx_tuples(setup_comb(1), 1), setup_idx_tuples(setup_comb(2), 1)};
        % Run kernel density estimate
        for j=1:2
            [f, xi] = ksdensity(feat_pair_mat_cell{j}, [ks_support_mesh_x(:), ks_support_mesh_y(:)]);
            if min(f) < caxis_common(1)
                caxis_common(1) = min(f);
            end
            if max(f) > caxis_common(2)
                caxis_common(2) = max(f);
            end
        end
        caxis_common_cell{setup_idx_tuples(setup_comb(1), 1), setup_idx_tuples(setup_comb(2), 1)} = caxis_common;
    end
end

% Draw the plots
% c_map = parula(length(plot_periods));
for i=1:size(setup_combs, 1) % Iterate over different combinations of statistic
    setup_comb = setup_combs(i, :);
    % do not include var and delta in the unification process
    if any(ismember(setup_comb, [5:5:20, 4:5:19]))
        continue; 
    end
    
    if ismember(2, setup_comb) || ismember(3, setup_comb)
        xaxis_opt = 'dist'; 
        xaxis_scale_opt = 10;
    else
        xaxis_opt = 'angle'; 
        xaxis_scale_opt = pi;
    end
    
    if setup_comb(2) == 27 || setup_comb(2) == 28
        yaxis_opt = 'attacked';
    else
        yaxis_opt = 'attacker';
    end
    
    ax_objs = cell(2, 1);
    for l=1:length(period_options)
        chosen_period = period_options(l);
        if ~ismember(chosen_period, plot_periods)
            continue;
        end
    
        feat_pair_mat_cell = cell(2, 1);
        for j=1:2 % Iterate over wildtype and mutant
            feat_pair_mat = [];
            for k=1:length(setup_comb)
                setup_tuple = setup_idx_tuples(setup_comb(k), :);
                feat_pair_mat(:, k) = featAll_cell{setup_tuple(1), setup_tuple(2), l, chosen_is_attacked_fly_idx, j};    
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
            
            ax_objs{j} = subplot(1, 2, j);
            surf(reshape(xi(:,1), ks_num_points, ks_num_points), reshape(xi(:,2), ks_num_points, ks_num_points), reshape(f, ks_num_points, ks_num_points), 'EdgeColor', 'none');
            colormap(teals); 
            view(2);
            
            axis tight;
            for k = 1:length(setup_comb)
                setup_tuple = setup_idx_tuples(setup_comb(k), :);
                if k == 1
%                     xlabel(strrep(strcat(feature_options{setup_tuple(1)}, '-', stat_options{setup_tuple(2)}), '_', '-'));
                    if strcmp(xaxis_opt, 'angle')
                        xlabel('Facing angle of attacking fly');
                    else
                        xlabel('Distance between two flies');
                    end
                    if contains(feature_options{setup_tuple(1)}, 'angle')
                        xticks((0:3).*(pi/3));
                        xticklabels({'0', '\pi/3', '2\pi/3', '\pi'});
                        xlim([0, pi]);
                    end
                else
%                     ylabel(strrep(strcat(feature_options{setup_tuple(1)}, '-', stat_options{setup_tuple(2)}), '_', '-'));
                    if strcmp(xaxis_opt, 'angle')
                        ylabel('Facing angle of attacked fly');
                    elseif strcmp(yaxis_opt, 'attacked')
                        ylabel('Facing angle of attacked fly');
                    else
                        ylabel('Facing angle of attacking fly');
                    end
                        
                    if contains(feature_options{setup_tuple(1)}, 'angle')
                        yticks((0:3).*(pi/3));
                        yticklabels({'0', '\pi/3', '2\pi/3', '\pi'});
                        ylim([0, pi]);
                    end
                end
            end
            
            quadrant_count = quadrant_count_cell{l, i, j};
            for m = 1:size(quadrant_count, 1)
                for n = 1:size(quadrant_count, 2)
                    text(m*xaxis_scale_opt/size(quadrant_count, 1) - 0.5*xaxis_scale_opt/size(quadrant_count, 1), n*pi/size(quadrant_count, 2) - 0.5*pi/size(quadrant_count, 1), 1, ...
                        {sprintf('  %d%%', round(100*quadrant_count(m, n)/sum(sum(quadrant_count)))), sprintf('N = %d', quadrant_count(m, n))});
                end
                if m ~= size(quadrant_count, 1)
                    line([0, xaxis_scale_opt], [m*pi/size(quadrant_count, 1), m*pi/size(quadrant_count, 1)], [1, 1], 'LineStyle', '--');
                    line([m*xaxis_scale_opt/size(quadrant_count, 1), m*xaxis_scale_opt/size(quadrant_count, 1)], [0, pi], [1, 1], 'LineStyle', '--');
                end
            end

            title(strrep(genotype_str{j}, '_', '-'));
        end

        caxis_common = caxis_common_cell{setup_idx_tuples(setup_comb(1), 1), setup_idx_tuples(setup_comb(2), 1)};
        for j=1:2
            caxis(ax_objs{j}, caxis_common);
            if j == 2
                original_pos = get(ax_objs{j}, 'Position');
                colorbar; 
                set(ax_objs{j}, 'Position', original_pos);
            end
        end
        % linkaxes([ax_objs{:}], 'xy');
        set(fig_obj, 'Position', get(0, 'Screensize'));
        suptitle(sprintf('Period %d', chosen_period));

        if chosen_is_attacked_fly
            fly_str = 'attacked';
        else
            fly_str = 'attacking';
        end
        saveas(double(fig_obj), sprintf('teal_heat_map-wt_vs_mutant-%s_vs_%s-period_%d-%s_fly.png', ...
            get(get(gca, 'ylabel'), 'String'), get(get(gca, 'xlabel'), 'String'), chosen_period, fly_str));
        saveas(double(fig_obj), sprintf('teal_heat_map-wt_vs_mutant-%s_vs_%s-period_%d-%s_fly.eps', ...
            get(get(gca, 'ylabel'), 'String'), get(get(gca, 'xlabel'), 'String'), chosen_period, fly_str));
    end         
end