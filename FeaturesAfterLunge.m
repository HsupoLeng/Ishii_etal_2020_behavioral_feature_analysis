function [feat_probs, feat_prob_edges, featAll] = FeaturesAfterLunge(flymat_name, common_path, feat_or_trk, feat_name, num_frames_vec, check_attacked_fly, remove_short_interval, remove_outliers, genotypes, selected_genotype, hist_stat)
    load(fullfile(common_path, flymat_name));
    featAll = struct('movie', '', 'fly', nan);

    if ~isempty(genotypes) % To look at pairing of specific genotypes, use, for example, [10,11] or [90, 1]; otherwise, [].
        all_genotypes = [flymatAll(:).genotype];
        flymatAll_mask_cell = arrayfun(@(type) all_genotypes == type, genotypes, 'UniformOutput', false);
        flymatAll_mask = any(vertcat(flymatAll_mask_cell{:}));
        flymatAll = flymatAll(flymatAll_mask);
    end
    
    featAll_sz = 0;
    for i = 1:length(flymatAll)
        % Allow mixed pairing. Use selected_genotype to indicate the
        % genotype that you want to examine in a pair
        if isempty(selected_genotype)
            selected_genotype = genotypes;
        end
        if ~ismember(flymatAll(i).genotype, selected_genotype) 
            continue; 
        end
        
        movie = flymatAll(i).movie{1};
        fly_feature = flymatAll(i).fly; % Take feature value from fly_feature
        another_fly_in_pair = fly_feature + (-1)^(mod(fly_feature, 2)+1);
        if check_attacked_fly % If we are interested in the feature of the attacked fly
            fly_time = another_fly_in_pair; 
        else
            fly_time = fly_feature; 
        end
        
        fly_time_idx = find(bitand(strcmp([flymatAll(:).movie], movie), [flymatAll(:).fly] == fly_time)); 
        LungeStarts = flymatAll(fly_time_idx).L_startsm;
        LungeEnds = flymatAll(fly_time_idx).L_endsm;
        
        entering_frame_30_min = min(flymatAll([fly_feature, another_fly_in_pair]).EnteringFrame);
        ending_frame_30_min = entering_frame_30_min + flymatAll(fly_feature).ThirtyMinFrame; 
        lunge_starts_mask = bitand(LungeStarts >= entering_frame_30_min, LungeStarts < ending_frame_30_min);
        lunge_ends_mask = bitand(LungeEnds >= entering_frame_30_min, LungeEnds < ending_frame_30_min);
        lunge_joint_mask = bitand(lunge_starts_mask, lunge_ends_mask);
        LungeStarts = LungeStarts(lunge_joint_mask);
        LungeEnds = LungeEnds(lunge_joint_mask);

        file_path = fullfile(common_path, movie(1:6), movie, movie);
        if ~feat_or_trk
            feat_mat_suffix = '-feat.mat';
            feat_struct_name = 'feat';
        else
            feat_mat_suffix = '-track.mat';
            feat_struct_name = 'trk';
        end
        feat_mat_name = fullfile(file_path, strcat(movie, feat_mat_suffix));

        try
            s_temp = load(feat_mat_name);
            feat = s_temp.(feat_struct_name);
        catch ME
            movie_parts = strsplit(movie, {'-', '_'});
            cwd_contents = dir(common_path);
            cwd_subfolders = {cwd_contents([cwd_contents(:).isdir]).name};
            target_subfolder = cwd_subfolders(contains(cwd_subfolders, movie_parts{1}));
            file_path = fullfile(common_path, target_subfolder{1}, movie, movie);
            feat_mat_name = fullfile(file_path, strcat(movie, feat_mat_suffix));
            s_temp = load(feat_mat_name);
            feat = s_temp.(feat_struct_name);
        end
        total_num_frame = size(feat.data, 2);
        
        if num_frames_vec(1) < 0
            LungeEnds(LungeStarts < abs(num_frames_vec(1))) = [];
            LungeStarts(LungeStarts < abs(num_frames_vec(1))) = [];
        end
        if num_frames_vec(2) > 0
            LungeStarts(abs(LungeEnds - total_num_frame) < num_frames_vec(2)) = [];
            LungeEnds(abs(LungeEnds - total_num_frame) < num_frames_vec(2)) = []; 
        end
        
        featAll_sz = featAll_sz + 1;
        featAll(featAll_sz).lunge_starts = LungeStarts';
        featAll(featAll_sz).lunge_ends = LungeEnds';
        featAll(featAll_sz).inter_lunge_interval = LungeStarts(2:end)' - LungeEnds(1:end-1)';
        
        if any(contains(feat_name, ' '))
            feat_name = cellfun(@(name) strrep(name, ' ', '_'), feat_name, 'UniformOutput', false);
        end
        if any(contains(feat.names, ' '))
            feat.names = cellfun(@(name) strrep(name, ' ', '_'), feat.names, 'UniformOutput', false);
        end
        
        mask_cell = cell(length(feat_name), 1);
        for j=1:length(feat_name)
            feature_str_components = strsplit(feat_name{j}, '_');
            feature_str_core = feature_str_components{1};
            mask_cell{j} = ~cellfun(@isempty, regexp(feat.names, strcat('^', feature_str_core)));
            cnt = 1;
            while sum(mask_cell{j}) > 1
                feature_str_core = strjoin({feature_str_core, feature_str_components{cnt + 1}}, '_');
                mask_cell{j} = ~cellfun(@isempty, regexp(feat.names, strcat('^', feature_str_core)));
                cnt = cnt + 1;
            end
        end
        
        featAll(featAll_sz).movie = movie; 
        featAll(featAll_sz).fly = fly_feature;
        featAll(featAll_sz).genotype = flymatAll(fly_feature).genotype; 

        for j=1:length(feat_name)
            if contains(feat_name{j}, 'mutual_other') 
                if mod(fly_time, 2)
                    fly_feature = fly_time + 1; % Take lunge times from fly_time
                else
                    fly_feature = fly_time - 1; 
                end
            else
                fly_feature = fly_time; 
            end
            if ~contains(feat_name{j}, 'lunge_interval')
                feat_vals = zeros(length(LungeStarts), diff(num_frames_vec) + 1);
                for k=1:length(LungeStarts)
                    feat_vals(k, :) = feat.data(fly_feature, LungeStarts(k)+num_frames_vec(1):LungeStarts(k)+num_frames_vec(2), mask_cell{j});
                end
            else
                feat_vals = cell(length(LungeStarts)-1, 1);
                for k=1:length(LungeStarts)-1
                    feat_vals{k} = feat.data(fly_feature, LungeEnds(k):LungeStarts(k+1)-1, mask_cell{j});
                end
            end
            featAll(featAll_sz).(feat_name{j}) = feat_vals;
        end
    end
    
    for i=1:length(featAll)
        for j=1:length(feat_name)
            if ~contains(feat_name{j}, 'lunge_interval')
                featAll(i).(strcat(feat_name{j}, '_mean')) = mean(featAll(i).(feat_name{j}), 2);
                featAll(i).(strcat(feat_name{j}, '_var')) = var(featAll(i).(feat_name{j}), 0, 2);
                featAll(i).(strcat(feat_name{j}, '_init')) = featAll(i).(feat_name{j})(:, 1);
                featAll(i).(strcat(feat_name{j}, '_end')) = featAll(i).(feat_name{j})(:, end);
                featAll(i).(strcat(feat_name{j}, '_delta')) = featAll(i).(feat_name{j})(:, end) - featAll(i).(feat_name{j})(:, 1);
                featAll(i).(strcat(feat_name{j}, '_max')) = max(featAll(i).(feat_name{j}), [], 2);
                featAll(i).(strcat(feat_name{j}, '_timepoint')) = featAll(i).(strcat(feat_name{j}, '_end'));
            else
                if isempty(featAll(i).(feat_name{j}))
                    continue; 
                else
                    featAll(i).(strcat(feat_name{j}, '_mean')) = cellfun(@(vec) mean(vec), featAll(i).(feat_name{j}));
                    featAll(i).(strcat(feat_name{j}, '_var')) = cellfun(@(vec) var(vec), featAll(i).(feat_name{j}));
                    feature_between_lunge_padded = cellfun(@(vec) [vec, nan(1, length(length(vec):diff(num_frames_vec)))], featAll(i).(feat_name{j}), 'UniformOutput', false);
                    featAll(i).(strcat(feat_name{j}, '_init')) = cellfun(@(vec) vec(1), feature_between_lunge_padded);
                    featAll(i).(strcat(feat_name{j}, '_end')) = cellfun(@(vec) vec(length(vec)), feature_between_lunge_padded);
                    featAll(i).(strcat(feat_name{j}, '_delta')) = cellfun(@(vec) vec(end)-vec(1), featAll(i).(feat_name{j}));
                    featAll(i).(strcat(feat_name{j}, '_max')) = cellfun(@max, featAll(i).(feat_name{j}));
                    if all(num_frames_vec < 0)
                        featAll(i).(strcat(feat_name{j}, '_timepoint')) = featAll(i).(strcat(feat_name{j}, '_init'));
                    else
                        featAll(i).(strcat(feat_name{j}, '_timepoint')) = cellfun(@(vec) vec(num_frames_vec(2)), feature_between_lunge_padded);
                    end
                end
            end
        end
    end
    
    % Remove sequences that have
    % more lunge happening in the duration designated by num_frames    
    featAll_fields = fieldnames(featAll);
    featAll_fields(contains(featAll_fields, 'lunge_interval')) = [];
    for i=1:length(featAll)
        if remove_short_interval
            if isempty(featAll(i).inter_lunge_interval)
                continue;
            end
            period_mask = false(size(featAll(i).lunge_starts));
            if num_frames_vec(1) < 0
                period_mask = [false; featAll(i).inter_lunge_interval < abs(num_frames_vec(1))];
            end
            if num_frames_vec(2) > 0
                period_mask = bitor(period_mask, [featAll(i).inter_lunge_interval < num_frames_vec(2); false]);
            end

            for k=1:length(featAll_fields)
                if contains(featAll_fields{k}, {'movie', 'fly', 'genotype', 'table'})
                    continue; 
                else
                    featAll(i).(featAll_fields{k})(period_mask, :) = [];
                end
            end
        end
    end
    
    % Remove outlier sequences using quartile method
    if all(num_frames_vec > 0)
        opt_str = 'post-lunge';
    elseif all(num_frames_vec < 0)
        opt_str = 'pre-lunge';
    else
        opt_str = 'peri-lunge';
    end
    if remove_outliers
        num_of_outliers = 0;
        for j=1:length(feat_name)
            if contains(feat_name{j}, {'angle', 'mutual'})
                continue;
            end
            curr_feature_all_means = vertcat(featAll(:).(strcat(feat_name{j}, '_mean'))); 
            quantiles = quantile(curr_feature_all_means, [0.005, 0.995]);
            outlier_mask = bitor(curr_feature_all_means < quantiles(1), curr_feature_all_means > quantiles(2));
            num_of_outliers = num_of_outliers + sum(outlier_mask);
            outlier_mask_by_fly_cell = mat2cell(outlier_mask, arrayfun(@(s) length(s.(strcat(feat_name{j}, '_mean'))), featAll));
            for i=1:length(outlier_mask_by_fly_cell)
                for k=1:length(featAll_fields)
                    if contains(featAll_fields{k}, {'movie', 'fly', 'genotype', 'table'})
                        continue; 
                    else
                        featAll(i).(featAll_fields{k})(outlier_mask_by_fly_cell{i}, :) = [];
                    end
                end
            end
        end
        fprintf('Removed %d outliers; %d %s sequences remaining\n', num_of_outliers, length(vertcat(featAll(:).(strcat(feat_name{j}, '_mean')))), opt_str);
    else
        fprintf('Did not remove outliers; %d %s sequences remaining\n', length(vertcat(featAll(:).(strcat(feat_name{j}, '_mean')))), opt_str);
    end
    
    if check_attacked_fly
        fly_str = 'attacked';
    else
        fly_str = 'attacking';
    end
    
    feat_probs = cell(length(feat_name), length(hist_stat));
    feat_prob_edges = cell(length(feat_name), length(hist_stat));
    for i=1:length(feat_name)
        if contains(feat_name{i}, 'lunge_interval')
            continue; 
        end
        for j=1:length(hist_stat)
            feat_stat = {featAll(:).(strcat(feat_name{i}, '_', hist_stat{j}))};
            feat_stat_zscore = zscore(vertcat(feat_stat{:}));
            if regexp(feat_name{i}, 'angle')
                bin_width = 0.1; 
            elseif regexp(feat_name{i}, 'dist')
                bin_width = 1;
            elseif regexp(feat_name{i}, 'vel')
                bin_width = 10; 
            end

            if regexp(hist_stat{j}, 'var')
                bin_width = bin_width^2; 
            end
            [feat_probs{i, j}, feat_prob_edges{i, j}] = histcounts(vertcat(feat_stat_zscore), 'BinWidth', bin_width, 'Normalization', 'probability');
        end
    end
end