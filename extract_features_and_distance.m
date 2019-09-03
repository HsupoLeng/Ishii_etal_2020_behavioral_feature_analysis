% All the options on features, statistics, time period that we use to
% collect the data
feature_options = {'dist_to_other', 'facing_angle', 'angle_between', 'vel', 'facing_angle_mutual_self', 'facing_angle_mutual_other'};
stat_options = {'mean', 'init', 'end', 'delta', 'var'};
period_options = [-60, -30, -10, -5, -3, -1, 10, 20, 30, 60];
is_attacked_fly = false;
fps = 60;

% Collect original features, feature statistics and compute the distance
% between wildtype and mutant feature statistic distributions
feat_probs_all_cell = cell(length(feature_options), length(stat_options), length(period_options), length(is_attacked_fly), 2);
feat_probs_dist_hellinger = zeros(length(feature_options), length(stat_options), length(period_options), length(is_attacked_fly));
feat_probs_dist_ks = zeros(length(feature_options), length(stat_options), length(period_options), length(is_attacked_fly));
featAll_cell = cell(length(feature_options), length(stat_options), length(period_options), length(is_attacked_fly), 2);
period_mask_all_cell = cell(length(period_options), length(is_attacked_fly), 2);

flymat_info_struct = struct();
% Same genotype pairing
flymat_info_struct(1).exp_folder = {'Z:\Kenichi\CSMH_SM_HeisenbergChamber_Apr-May2019', 'Z:\Kenichi\CG3385_KO_CS-bc-11_HeisenbergChamber_Apr-May2019'};
flymat_info_struct(1).flymat = {'FLYMAT_CSMH_SM_HeisenbergChamber_Apr-May2019', 'FLYMAT_CG3385_KO_CS-bc-11_HeisenbergChamber_Apr2019'};
flymat_info_struct(1).genotypes = {[90, 91], [10, 11]};
flymat_info_struct(1).selected_genotype = {[], []};
flymat_info_struct(1).genotype_str = {'WT_SM', 'CG3385_KO_CS-bc-11_GM'};

% SH WT vs. mutant
flymat_info_struct(2).exp_folder = {'Z:\Kenichi\CG3385_KO_CS-bc-11_GM_vsWT_SM_HeisenbergChamber_May2019'};
flymat_info_struct(2).flymat = {'FLYMAT_CG3385_KO_CS-bc-11_GM_vsWT_SM_HeisenbergChamber_May2019'};
flymat_info_struct(2).genotypes = {[90, 92, 10, 12]};
flymat_info_struct(2).selected_genotype = {[90, 92], [10, 12]};
flymat_info_struct(2).genotype_str = {'WT_SM', 'CG3385_KO_CS-bc-11_GM'};

% GH WT vs. SH WT
flymat_info_struct(3).exp_folder = {'Z:\Kenichi\CSMH_SM_vsWT_GM_HeisenbergChamber_May2019'};
flymat_info_struct(3).flymat = {'FLYMAT_CSMH_SM_vsWT_GM_HeisenbergChamber_May2019'};
flymat_info_struct(3).genotypes = {[90, 91]};
flymat_info_struct(3).selected_genotype = {[90], [91]};
flymat_info_struct(3).genotype_str = {'WT_GM', 'WT_SM'};

% GH WT vs. mutant
flymat_info_struct(4).exp_folder = {'Z:\Kenichi\CG3385_KO_CS-bc-11_GM_vsWT_GM_HeisenbergChamber_May2019'};
flymat_info_struct(4).flymat = {'FLYMAT_CG3385_KO_CS-bc-11_GM_vsWT_GM_HeisenbergChamber_May2019'};
flymat_info_struct(4).genotypes = {[1, 90]};
flymat_info_struct(4).selected_genotype = {[90], [1]};
flymat_info_struct(4).genotype_str = {'WT_GM', 'CG3385_KO_CS-bc-11_GM'};

flymat_sel = 1;
for k=1:length(period_options)
    for l=1:length(is_attacked_fly)
        fprintf('Now at set-up %d/%d\n', (k-1)*length(is_attacked_fly)+l, length(period_options)*length(is_attacked_fly));
          
        [feat_probs_wt_cell, feat_prob_edges_wt_cell, featAll_wt] = FeaturesAfterLunge(flymat_info_struct(flymat_sel).flymat{1}, flymat_info_struct(flymat_sel).exp_folder{1}, feature_options, period_options(k), is_attacked_fly(l), flymat_info_struct(flymat_sel).genotypes{1}, flymat_info_struct(flymat_sel).selected_genotype{1}, stat_options);
        if length(flymat_info_struct(flymat_sel).flymat) == 1
            [feat_probs_mu_cell, feat_prob_edges_mu_cell, featAll_mu] = FeaturesAfterLunge(flymat_info_struct(flymat_sel).flymat{1}, flymat_info_struct(flymat_sel).exp_folder{1}, feature_options, period_options(k), is_attacked_fly(l), flymat_info_struct(flymat_sel).genotypes{1}, flymat_info_struct(flymat_sel).selected_genotype{2}, stat_options);
        else
            [feat_probs_mu_cell, feat_prob_edges_mu_cell, featAll_mu] = FeaturesAfterLunge(flymat_info_struct(flymat_sel).flymat{2}, flymat_info_struct(flymat_sel).exp_folder{2}, feature_options, period_options(k), is_attacked_fly(l), flymat_info_struct(flymat_sel).genotypes{2}, flymat_info_struct(flymat_sel).selected_genotype{2}, stat_options);
        end
        % Store the features
        for i=1:length(feature_options)
            for j=1:length(stat_options)
                field_name = strcat(feature_options{i}, '_', stat_options{j});
                featAll_cell{i, j, k, l, 1} = vertcat(featAll_wt(:).(field_name));
                featAll_cell{i, j, k, l, 2} = vertcat(featAll_mu(:).(field_name));
            end
        end
        
%         % Generate period_mask
%         period_mask_all_cell{k, l, 1} = logical(vertcat(featAll_wt(:).period_mask));
%         period_mask_all_cell{k, l, 2} = logical(vertcat(featAll_mu(:).period_mask));
        
        % Compute statistical distance between wildtype and mutant feature distributions
        for i=1:length(feature_options)
            for j=1:length(stat_options)
                feat_probs_wt = feat_probs_wt_cell{i, j};
                feat_probs_mu = feat_probs_mu_cell{i, j};
                feat_prob_edges_wt = feat_prob_edges_wt_cell{i, j};
                feat_prob_edges_mu = feat_prob_edges_mu_cell{i, j};

                if min(feat_prob_edges_wt) ~= min(feat_prob_edges_mu)
                    if min(feat_prob_edges_wt) < min(feat_prob_edges_mu)
                        feat_probs_mu = padarray(feat_probs_mu, round([0, (min(feat_prob_edges_mu) - min(feat_prob_edges_wt))/unique(round(diff(feat_prob_edges_mu), 2))]), 'pre');
                    else
                        feat_probs_wt = padarray(feat_probs_wt, round([0, (min(feat_prob_edges_wt) - min(feat_prob_edges_mu))/unique(round(diff(feat_prob_edges_mu), 2))]), 'pre');
                    end
                end
                if max(feat_prob_edges_wt) ~= max(feat_prob_edges_mu)
                    if max(feat_prob_edges_wt) < max(feat_prob_edges_mu)
                        feat_probs_wt = padarray(feat_probs_wt, round([0, (max(feat_prob_edges_mu) - max(feat_prob_edges_wt))/unique(round(diff(feat_prob_edges_mu), 2))]), 'pre');
                    else
                        feat_probs_mu = padarray(feat_probs_mu, round([0, (max(feat_prob_edges_wt) - max(feat_prob_edges_mu))/unique(round(diff(feat_prob_edges_mu), 2))]), 'pre');
                    end
                end
                feat_probs_wt_cell{i, j} = feat_probs_wt; 
                feat_probs_mu_cell{i, j} = feat_probs_mu;
            end
        end

        feat_probs_all_cell(:, :, k, l, 1) = feat_probs_wt_cell;
        feat_probs_all_cell(:, :, k, l, 2) = feat_probs_mu_cell; 
        feat_probs_dist_hellinger(:, :, k, l) = cellfun(@(feat_probs_wt, feat_probs_mu) sqrt(1 - (sum(sqrt(feat_probs_wt.*feat_probs_mu)))), feat_probs_wt_cell, feat_probs_mu_cell); % Compute Hellinger distance
        feat_probs_dist_ks(:, :, k, l) = cellfun(@(feat_probs_wt, feat_probs_mu) max(abs(feat_probs_wt - feat_probs_mu)), feat_probs_wt_cell, feat_probs_mu_cell); % Compute Kolmogorov-Smirnov distance
    
        % Plot histogram for inter-lunge interval
        % The histogram can be different for different periods due to
        % different outliers removed. Currently we draw figures based on data
        % of the attacked flym with outliers removed using 10-frame
        % sequences
        if k == 1 && l == 1
            ili_wt = vertcat(featAll_wt(:).inter_lunge_interval)./fps; 
            ili_mu = vertcat(featAll_mu(:).inter_lunge_interval)./fps;
            figure();
            histogram(ili_wt, 'BinWidth', 30);
            title('Inter-lunge interval histogram for wild-type flies');
            xlabel('Inter-lunge interval (in second)');
            ylabel('Number of occurrences');
            saveas(double(gcf), 'histogram_inter_lunge_interval_wt_all.eps');
            figure();
            histogram(ili_wt(ili_wt<10), 'BinWidth', 0.5);
            saveas(double(gcf), 'histogram_inter_lunge_interval_wt_zoom_in.eps');
            figure();
            histogram(ili_mu, 'BinWidth', 30);
            title('Inter-lunge interval histogram for mutant flies');
            xlabel('Inter-lunge interval (in second)');
            ylabel('Number of occurrences');
            saveas(double(gcf), 'histogram_inter_lunge_interval_mu_all.eps');
            figure();
            histogram(ili_mu(ili_mu<10), 'BinWidth', 0.5);
            saveas(double(gcf), 'histogram_inter_lunge_interval_mu_zoom_in.eps');
        end
    end
end

genotype_str = flymat_info_struct(flymat_sel).genotype_str; 
% save the result
save(strcat('featAll_and_feat_probs_dist_remove_outliers-', strjoin(flymat_info_struct(flymat_sel).flymat, '_'), '.mat'), 'feature_options', 'stat_options', 'period_options', 'is_attacked_fly', 'featAll_cell', 'feat_probs_all_cell', 'feat_probs_dist_hellinger', 'feat_probs_dist_ks', 'period_mask_all_cell', 'genotype_str');