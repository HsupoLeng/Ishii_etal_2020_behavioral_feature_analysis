function [info, cluster_labels_cell_by_bout, cluster_bout_length_cell, cluster_labels_changepoint_cell, labels_all] = ...
    apply_behav_label(base_path, flymat_file, genotype_select, genotype_name, remove_fragment)
% example input 1:
% base_path = 'D:\xubo\behavioral_movies\CG3385KO_CS-bc-11_Df-Exel6082_CS-bc-6_GM_Aug2019'; 
% flymat_file = 'FLYMAT_CG3385KO_CS-bc-11_Df-Exel6082_CS-bc-6_GM_Aug2019.mat'' 
% genotype_select = [1, 4];
% genotype_name = {'wildtype GH', 'delta-nvy GH'};

% example input 2:
% base_path = 'D:\xubo\behavioral_movies\CG3385_KO_CS-bc-11_SM_Mar2017';
% flymat_file = 'FLYMAT_CG3385_KO_CS-bc-11_SM_Mar2017.mat';
% genotype_select = [99]; 
% genotype_name = {'wildtype SH'};

    fps = 60; 
    speed_threshold_low = 1; % Robie et.al., (JEB, 2010)
    speed_threshold_high = 2.5; % 5mm/s would exclude smaller walks
    dist_threshold_low = 5.0; 
    dist_threshold_high = 7.5; 
    angle_threshold_low = pi/4;
    angle_threshold_high = pi/3;
    num_of_label = 6; 
    behav_name = {'stopping', 'orienting', 'non-orienting', 'lunging', 'wing extension', 'fragment'};
    jaaba_smoothed = true; 
    max_gap = 3; 
    min_bout = 6;
    if ~exist('remove_fragment', 'var')
        remove_fragment = false; 
    end

    cluster_labels_cell_by_bout = cell(size(genotype_select)); 
    cluster_bout_length_cell = cell(size(genotype_select)); 
    cluster_labels_changepoint_cell = cell(size(genotype_select)); 

    speed_cell = cell(size(genotype_select)); 
    dist_cell = cell(size(genotype_select)); 
    fa_angle_cell = cell(size(genotype_select)); 
    
    if ischar(base_path)
        load(fullfile(base_path, flymat_file)); 
    else
        assert(iscell(base_path)); 
        assert(iscell(flymat_file)); 
        assert(length(base_path) == length(genotype_select)); 
    end
    for g=1:length(genotype_select)
        if ~ischar(base_path)
            load(fullfile(base_path{g}, flymat_file{g})); 
        end
        fly_idx = ismember([flymatAll.genotype], genotype_select{g}); 
        flymat_temp = flymatAll(fly_idx);
        movie_name = unique([flymatAll(fly_idx).movie]); 

        movie_date = cellfun(@(s) strsplit(s, '_'), movie_name, 'UniformOutput', false); 
        movie_date = cellfun(@(c) c{1}, movie_date, 'UniformOutput', false); 
        for m=1:length(movie_date)
            if contains(movie_date{m}, '170718')
                movie_date{m} = '170718-1';
            end
            if strcmp(flymat_file, 'FLYMAT_CSMH_SM_vsWT_GM_HeisenbergChamber_May2019.mat')  
                movie_date{m} = strcat(movie_date{m}, '-2'); 
            end
            if strcmp(flymat_file, 'FLYMAT_Tdc2-GAL4_nvy-LexA_Tub-FRT-GAL80_LexAop-FLPL_UAS-CsChrim_SM_vsCSmale_LED5minOffOnOff_Feb2019')
                movie_date{m} = strcat(movie_date{m}, '_opto'); 
            end
            if contains(movie_date{m}, '190517') && strcmp(flymat_file, 'FLYMAT_Tdc2-GAL4_UAS-IR-inaF_CG33172_GstE12_GM_May2019')
                movie_date{m} = strcat(movie_date{m}, '-2'); 
            end
        end

        if ischar(base_path)
            feat_path = cellfun(@(d, n) fullfile(base_path, d, n, n, strcat(n, '_JAABA')), movie_date, movie_name, 'UniformOutput', false); 
        else
            feat_path = cellfun(@(d, n) fullfile(base_path{g}, d, n, n, strcat(n, '_JAABA')), movie_date, movie_name, 'UniformOutput', false); ; 
        end
        x_mm = [];
        y_mm = [];
        speed = [];
        speed_ang = []; 
        dist = [];
        fa_angle = []; 

        for m=1:length(feat_path)
            fly_idx_one_mov = find_fly_in_flymat(flymat_temp, movie_name{m}, [], true); 
            fly_rel_idx = [flymat_temp(fly_idx_one_mov).fly]; 

            cd(feat_path{m}); 
            load('trx.mat'); 

            x_mm_temp = {trx(fly_rel_idx).x_mm};
            y_mm_temp = {trx(fly_rel_idx).y_mm};
            x_mm = horzcat(x_mm, x_mm_temp); 
            y_mm = horzcat(y_mm, y_mm_temp);  
            clear trx;


            cd('.\perframe'); 
            fa_angle_temp = load('facing_angle.mat', 'data');
            log_speed_temp = load('log_vel.mat', 'data');
            dist_to_other_temp = load('dist_to_other.mat', 'data'); 
            speed_temp = cell(size(fly_rel_idx)); 
            speed_ang_temp = cell(size(fly_rel_idx)); 
            dist_temp = cell(size(fly_rel_idx)); 
            for i = 1:length(speed_temp)
                speed_temp{i} = 10.^log_speed_temp.data{fly_rel_idx(i)}; 
                dist_temp{i} = dist_to_other_temp.data{fly_rel_idx(i)}; 
%                 speed_temp{i} = 10.^log_speed_temp.data{i}; %to-do:
%                 determine if this was a bug for previous data
%                 dist_temp{i} = dist_to_other_temp.data{i}; 
    %             speed_temp{i} = sqrt(diff(x_mm_temp{i}).^2 + diff(y_mm_temp{i}).^2)./dt; 
                speed_ang_temp{i} = atan(-diff(y_mm_temp{i})./diff(x_mm_temp{i})); 
    %             if rem(i, 2) == 0 
    %                 d = (x_mm_temp{i-1} - x_mm_temp{i}).^2 + (y_mm_temp{i-1} - y_mm_temp{i}).^2; 
    %                 dist_temp{i-1} = ; 
    %                 dist_temp{i} = d; 
    %             end
            end
            speed = horzcat(speed, speed_temp); 
            speed_ang = horzcat(speed_ang, speed_ang_temp); 
            dist= horzcat(dist, dist_temp); 
            fa_angle = horzcat(fa_angle, fa_angle_temp.data(fly_rel_idx)); 
        end
        speed_cell{g} = speed; 
        speed_ang_cell{g} = speed_ang; 
        dist_cell{g} = dist; 
        fa_angle_cell{g} = fa_angle; 
    end
    clear x_mm x_mm_temp y_mm y_mm_temp speed speed_temp log_speed_temp speed_ang speed_ang_temp dist dist_to_other_temp dist_temp fa_angle fa_angle_temp;
    
    %% Determine non-interacting, stopping, non-orienting and orienting
    labels_all = cell(size(genotype_select)); 
    for g=1:length(genotype_select)
        if ~ischar(base_path)
            load(fullfile(base_path{g}, flymat_file{g})); 
        end
        fly_idx = ismember([flymatAll.genotype], genotype_select{g}); 
        labels = cell(size(speed_cell{g})); 
        for l=1:numel(labels)
            labels{l} = zeros(size(speed_cell{g}{l})); 
        end

        % stopping as base state, using speed as the only criterion
        stopping = cellfun(@(s) ~schmitt_trigger(s, speed_threshold_low, speed_threshold_high), speed_cell{g}, 'UniformOutput', false); 
        labels = cellfun(@(la, ni) la + ni, labels, stopping, 'UniformOutput', false); 
        clear stopping;

        orienting = cellfun(@(ang) 2.*(~schmitt_trigger(abs(ang), angle_threshold_low, angle_threshold_high)), fa_angle_cell{g}, 'UniformOutput', false); 
        orienting = cellfun(@(o, d) o.*(~schmitt_trigger(d, dist_threshold_low, dist_threshold_high)), orienting, dist_cell{g}, 'UniformOutput', false);
        labels = cellfun(@(la, o) la + o.*~la, labels, orienting, 'UniformOutput', false); 
        clear orienting;

        non_orienting = cellfun(@(la) 3.*~la, labels, 'UniformOutput', false); 
        labels = cellfun(@(la, no) la + no.*~la, labels, non_orienting, 'UniformOutput', false);
        clear non_orienting;


        assert(all(cellfun(@(la) sum(~la(:)) == 0, labels))); 
        assert(all(cellfun(@(la) sum(la(:)>3) == 0, labels))); 

        %% Merge labels with other labels from JAABA, e.g. lunge and wing extension
        behav_shorthand = {'L', 'WE'};
        behav_label = [4, 5]; 
        fly_idx = find(fly_idx);
        for fl=1:length(fly_idx)
            labels_jaaba = zeros(size(labels{fl})); 
            for b=1:length(behav_shorthand)
                if isfield(flymatAll, strcat(behav_shorthand{b}, '_binary'))
                    if jaaba_smoothed
                        behav_binary = flymatAll(fly_idx(fl)).(strcat(behav_shorthand{b}, '_binary'));
                        if contains(flymat_file, 'Heisenberg')
                            behav_binary = behav_binary(flymatAll(fly_idx(fl)).EnteringFrame:end); 
                            if length(behav_binary) ~= length(labels_jaaba)
                                if length(behav_binary) < length(labels_jaaba)
                                    behav_binary = [zeros(1, length(labels_jaaba) - length(behav_binary)), behav_binary]; 
                                else
                                    behav_binary = behav_binary((end-length(labels_jaaba)+1):end); 
                                end
                            end
                        end
                        if length(behav_binary) ~= length(labels_jaaba)
                            behav_binary = behav_binary(1:length(labels_jaaba));
                        end
                        labels_jaaba = labels_jaaba + behav_binary.*~labels_jaaba.*behav_label(b);
                    else
                        jaaba_bout_temp = arrayfun(@(s, e) s:e-1, flymatAll(fly_idx(fl)).(strcat(behav_shorthand{b}, '_start')), flymatAll(fly_idx(fl)).(strcat(behav_shorthand{b}, '_end')), 'UniformOutput', false);
                        if ~isempty(jaaba_bout_temp)
                            jaaba_binary_temp = accumarray(horzcat(jaaba_bout_temp{:})', 1, [length(flymatAll(fly_idx(fl)).(strcat(behav_shorthand{b}, '_binary'))), 1]); 
                        else
                            jaaba_binary_temp = zeros(size(flymatAll(fly_idx(fl)).(strcat(behav_shorthand{b}, '_binary'))))'; 
                        end
                        labels_jaaba = labels_jaaba + jaaba_binary_temp'.*~labels_jaaba.*behav_label(b);
                    end
                end
            end
            labels{fl} = labels{fl}.*~labels_jaaba + labels_jaaba; 
        end
        assert(all(cellfun(@(la) sum(~la(:)) == 0, labels))); 
        assert(all(cellfun(@(la) sum(la(:)>5) == 0, labels)));  

        %% Fill gap for orienting, non-orienting and stopping
        % Then label short fragment bouts as NUM_OF_LABEL
        bin_labels_changepoint_cell = cell(size(labels)); 
        bin_bout_length_cell = cell(size(labels)); 
        for b=1:3
            labels_bin = cellfun(@(l) l == b, labels, 'UniformOutput', false); 
            [bin_label_cell, bin_bout_length_cell, bin_labels_changepoint_cell] = ...
                convert_cluster_labels_to_bout(labels_bin, 0); 
            t0s = cellfun(@(p, la) p(la), bin_labels_changepoint_cell, bin_label_cell, 'UniformOutput', false); 
            bin_bout_length_cell  = cellfun(@(l, la) l(la), bin_bout_length_cell, bin_label_cell, 'UniformOutput', false);
            t1s = cellfun(@(t, l) t + l, t0s, bin_bout_length_cell, 'UniformOutput', false);
            allScores = struct();
            allScores.t0s = t0s; 
            allScores.t1s = t1s; 
            allScores.postprocessed = labels; % Not "post-processed" as it was in JAABA's case
            allScores = smoothing_xubo(allScores, 0, max_gap, min_bout);

            for f=1:length(labels)
                labels_one_fly = labels{f}; 
                old_mask = labels_bin{f}; 
                new_mask = logical(allScores.binary{f}); 
                labels_one_fly(new_mask) = b; 
                labels_one_fly(bitand(old_mask, ~new_mask)) = num_of_label; 
                labels{f} = labels_one_fly; 
            end
        end
        
        if remove_fragment
            for f=1:length(labels)
                labels_one_fly = labels{f}; 
                labels{f} = labels_one_fly(labels_one_fly ~= num_of_label); 
            end
        end
        
        
        %% Convert frames labels into bouts
        [cluster_labels_cell_by_bout{g}, cluster_bout_length_cell{g}, cluster_labels_changepoint_cell{g}] = ...
            convert_cluster_labels_to_bout(labels, 0); 
        
        labels_all{g} = labels; 
    end
    
    %% Compile label bout information
    info = struct();
    for g=1:length(genotype_select)
        if ~ischar(base_path)
            load(fullfile(base_path{g}, flymat_file{g})); 
        end
        info(g).genotype = genotype_name{g}; 
        fly_idx = find(ismember([flymatAll.genotype], genotype_select{g})); 
        flies = struct(); 
        for f=1:length(cluster_labels_cell_by_bout{g})
            flies(f).movie = flymatAll(fly_idx(f)).movie; 
            flies(f).fly = flymatAll(fly_idx(f)).fly; 

            flies(f).bouts = table(cluster_labels_changepoint_cell{g}{f}, ...
                cluster_labels_cell_by_bout{g}{f}, ...
                cluster_bout_length_cell{g}{f}, ...
                'VariableNames', {'startFrame', 'label', 'length'}); 
        end
        info(g).flies = flies; 
    end
end
