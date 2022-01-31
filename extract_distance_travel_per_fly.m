base_path = 'Z:\Kenichi';
output_path = 'D:\xubo\nervy_w_Kenichi\data-012222';

curr_movie_name = ''; 
remove_fragment = false; % has to be false for behavior label and frame-by-frame data to match.  
% non_orient_only = false; 
fps = 60; 
same_genotype_pairing = true;  


% ===== Fig. 1 =====
% exp_list = {'Elav-GAL4_UAS-IR-nvy(VDRC100373,JF03349)_May2018', ...
%     'CG3385KO_CS-bc-11_Df-Exel6082_CS-bc-6_GM_Aug2019', ...
%     'CG3385_KO_CS-bc-11_SM_Mar2017', ...
%     'CG3385KO(CS-bc-11)_Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_GM_Jul2017', ...
%     'Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_SM_Aug2017'};
% genotype_select = {[1, 2, 3], [1, 4], [99, 3], [2, 3, 1], [21, 31, 11]};
% genotype_name = {{'Elav-GAL4, 260b', '+, UAS-IR nvy', 'Elav-GAL4, UAS-IR nvy'}, {'wildtype GH', 'delta-nvy GH'}, ...
%     {'wildtype SH', 'delta-nvy SH'}, {'delta-nvy, elav-GAL4>GFP', 'delta-nvy, +>nvy', 'delta-nvy, elav-GAL4>nvy'}, ...
%     {'elav-GAL4>GFP', '+>nvy', 'elav-GAL4>nvy'}};
% time_period = {}; 
% movie_length_minute = 30;  

% ===== Fig. 3E =====
% exp_list = {'Tdc2-GAL4_nvy-LexA_Tub-FRT-GAL80-FRTorTub-FRT-stop-FRT-GAL80_LexAop-FLPL_UAS-Kir_GM_Nov2018', ...
%     'Tdc2-GAL4_nvy-LexA_Tub-FRT-GAL80-FRTorTub-FRT-stop-FRT-GAL80_LexAop-FLPL_UAS-Kir_SM_Nov2018'};
% genotype_select = {[23, 43, 63], [11, 21, 41, 61]};
% genotype_name = {{'Tdc2-GAL4, nvy-LexA, Tub>GAL80>, UAS-Kir2.1', 'Tdc2-GAL4, nvy-LexA, Tub>GAL80>, LexA-FLPL, UAS-GFP', 'Tdc2-GAL4, nvy-LexA, Tub>GAL80>, LexA-FLPL, UAS-Kir2.1'}, ...
%     {'Tdc2-GAL4, nvy-LexA, UAS-Kir2.1', 'Tdc2-GAL4, nvy-LexA, Tub>GAL80>, UAS-Kir2.1', 'Tdc2-GAL4, nvy-LexA, Tub>GAL80>, LexA-FLPL, UAS-GFP', 'Tdc2-GAL4, nvy-LexA, Tub>GAL80>, LexA-FLPL, UAS-Kir2.1'}};
% time_period = {}; 
% movie_length_minute = 30; 

% ===== Fig. 3G =====
% exp_list = {'Tdc2-GAL4_nvy-LexA_Tub-FRT-GAL80-FRTorTub-FRT-stop-FRT-GAL80_LexAop-FLPL_UAS-Kir_GM_Nov2018', ...
%     'Tdc2-GAL4_nvy-LexA_Tub-FRT-GAL80-FRTorTub-FRT-stop-FRT-GAL80_LexAop-FLPL_UAS-Kir_SM_Nov2018'};
% genotype_select = {[33, 53, 73], [31, 51, 71]};
% genotype_name = {{'Tdc2-GAL4, nvy-LexA, Tub>stop>GAL80, UAS-Kir2.1', 'Tdc2-GAL4, nvy-LexA, Tub>stop>GAL80, LexA-FLPL, UAS-GFP', 'Tdc2-GAL4, nvy-LexA, Tub>stop>GAL80, LexA-FLPL, UAS-Kir2.1'}, ...
%     {'Tdc2-GAL4, nvy-LexA, Tub>stop>GAL80, UAS-Kir2.1', 'Tdc2-GAL4, nvy-LexA, Tub>stop>GAL80, LexA-FLPL, UAS-GFP', 'Tdc2-GAL4, nvy-LexA, Tub>stop>GAL80, LexA-FLPL, UAS-Kir2.1'}};
% time_period = {}; 
% movie_length_minute = 30; 

% ===== Fig. S2C =====
% exp_list = {'CG3385KO(CS-bc-11)_elav-gal4,UAS-hMTGs(CS-bc-6)_Aug2017'};
% genotype_select = {[1, 4, 5, 6, 7]};
% genotype_name = {{'delta-nvy, elav-GAL4, UAS-GFP', 'delta-nvy, +, UAS-hMTG8', 'delta-nvy, elav-GAL4, UAS-hMTG8', 'delta-nvy, +, UAS-hMTG16', 'delta-nvy, elav-GAL4, UAS-hMTG16'}}; 
% time_period = {}; 
% movie_length_minute = 30; 

% ===== Fig. S2F =====
% exp_list = {'CG3385KO(CS-bc-11)_elav-gal4,UAS-nvy(CS-bc-6)_Mar,Jun2017'}; 
% genotype_select = {[1, 2, 3, 4, 5, 6, 8]}; 
% genotype_name = {{'delta-nvy, elav-GAL4, +', 'delta-nvy, elav-GAL4, GFP', 'delta-nvy, elav-GAL4, UAS-nvy', 'delta-nvy, elav-GAL4, UAS-nvy-delta-NHR1', 'delta-nvy, elav-GAL4, UAS-nvy-delta-NHR2', 'delta-nvy, elav-GAL4, UAS-nvy-delta-NHR3', 'delta-nvy, elav-GAL4, UAS-nvy-delta-NHR4'}}; 
% time_period = {}; 
% movie_length_minute = 30; 

% ===== Fig. 4 and S9 =====
exp_list = {'Tdc2-GAL4_UAS-IR-Yp2_pie_DNApole255_GM_Apr2019', ...
    'Tdc2-GAL4_UAS-IR-CG1552_GILT1_Cpsf73_Kua_GM_Apr2019', ...
    'Tdc2-GAL4_UAS-IR-CG4302_Cpsf100_CG3176_Dop1R2_GM_Apr2019', ...
    'Tdc2-GAL4_UAS-IR-CG8569_CG4038_CG14377_CG4554_CG3168_GM_Apr2019', ... 
    'Tdc2-GAL4_UAS-IR-inaF_CG33172_GstE12_GM_May2019', ...
    'Tdc2-GAL4_UAS-IR-CG12926_GstT2_GM_Sep2019', ...
    'Tdc2-GAL4_UAS-IR-IntS8_AstA_alrm_GM_Aug2019', ...
    'Tdc2-GAL4_UAS-IR-Irk2_Cyp9f2_CG5808_GM_Aug2019', ...
    'CG3385KO_Tdc2-GAL4_UAS-IR-Tango2_CG32113_GM_Jun2019', ... 
    'CG3385KO_Tdc2-GAL4_UAS-IR-CG7692_ZC3H3_GM_Jul2019', ...
    'CG3385KO_Tdc2-GAL4_UAS-IR-bmm_GM_Sep2019', ...
    'CG3385KO_Tdc2-GAL4_UAS-IR-CG8833_CG3294_CG18273_CG16771_GM_Oct2019', ...
    'CG3385KO_Tdc2-GAL4_UAS-IR-CG3587_CG18273_CG3071_GM_Oct2019'
    };
genotype_select = {[1:7], [1:9], [1:9], [1:11], [1:7], [1:5], [1:7], [1:7], [1:6], ...
    [1:6], [1:3], [1:8], [1:7]};
time_period = {}; 
movie_length_minute = 30; 

% ===== Optogenetic activation of nvy+/Tdc2+ =====
% exp_list = {'Tdc2-GAL4_nvy-LexA_Tub-FRT-GAL80_LexAop-FLPL_UAS-CsChrim_SM_vsCSmale_LED5minOffOnOff_Feb2019'}; 
% genotype_select = {[11, 21, 31]};
% genotype_name = {{'nvy-LexA, Tdc2-GAL4, Tub>GAL80>, +, UAS-CsChrimson', 'nvy-LexA, Tdc2-GAL4, Tub>GAL80>, LexAop2-FLPL, UAS-GFP', 'nvy-LexA, Tdc2-GAL4, Tub>GAL80>, LexAOp2-FLPL, UAS-CsChrimson'}}; 
% time_period = {'pre-opto', 'opto', 'post-opto'};
% time_period_length = 5*60*fps; 
% time_period_start = [0:length(time_period)-1]' .* time_period_length + 1;
% time_period_end = time_period_start + time_period_length - 1;
% movie_length_minute = 15; 

if same_genotype_pairing
    suffix = "per_pair";
else
    suffix = "per_fly"; 
end
dist_travel_struct = struct('exp', '', sprintf('dist_%s', suffix), [], sprintf('dist_%dminute_%s', round(movie_length_minute/3), suffix), [], sprintf('dist_minute_%s', suffix), [], sprintf('lungePerMeter_%s', suffix), nan, sprintf('speed_in_non_orient_%s', suffix), []); 
dist_travel_struct(1) = []; 

for e=1:length(exp_list)
    exp_path = fullfile(base_path, exp_list{e}); 
    cd(exp_path);
    try
        flymat_filename = strcat('FLYMAT_', exp_list{e}); 
        load(flymat_filename); 
    catch 
        flymat_file = dir('FLYMAT_*.mat');
        if strcmp(exp_path, 'Z:\Kenichi\Tdc2-GAL4_UAS-IR-CG12926_GstT2_GM_Sep2019')
            load('t.mat'); 
            flymat_filename = 't.mat'; 
        else
            load(flymat_file.name);
            flymat_filename = flymat_file.name; 
        end
    end
    
    [~, ~, ~, ~, labels_all] = ...
        apply_behav_label(exp_path, flymat_filename, {unique([flymatAll(:).genotype])}, {'mixed'}, remove_fragment);
    
    labels = labels_all{:}; 
    dist = nan(size(flymatAll)); 
    dist_onethird_portion = nan([length(flymatAll), 3]); 
    dist_minute = nan([length(flymatAll), movie_length_minute]); 
    if isempty(time_period)
        dist_in_non_orient = nan(size(flymatAll));
        time_in_non_orient = nan(size(flymatAll)); 
    else
        dist_in_non_orient = nan([length(flymatAll), length(time_period)]);
        time_in_non_orient = nan([length(flymatAll), length(time_period)]); 
    end
    
    for f=1:length(flymatAll)
        fly_id = flymatAll(f).fly; 
        movie_name = flymatAll(f).movie{1}; 
        movie_name_elem = strsplit(movie_name, '_'); 
        day_folder = movie_name_elem{1}; 
        if strcmp(exp_list{e}, 'CG3385KO(CS-bc-11)_Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_GM_Jul2017') && contains(movie_name, '170718')
            day_folder = strcat(day_folder, '-1'); 
        end
        if strcmp(exp_list{e}, 'Tdc2-GAL4_UAS-IR-inaF_CG33172_GstE12_GM_May2019') && contains(movie_name, '190517')
            day_folder = strcat(day_folder, '-2'); 
        end
        if strcmp(exp_list{e}, 'Tdc2-GAL4_nvy-LexA_Tub-FRT-GAL80_LexAop-FLPL_UAS-CsChrim_SM_vsCSmale_LED5minOffOnOff_Feb2019')
            day_folder = strcat(day_folder, '_opto'); 
        end
        if ~strcmp(movie_name, curr_movie_name)
            cd(fullfile(exp_path, day_folder, movie_name, movie_name, strcat(movie_name, '_JAABA')));       
            load('trx.mat');
            curr_movie_name = movie_name; 
        end
        
        tracked_mask = bitand(~isnan(trx(fly_id).x_mm), ~isnan(trx(fly_id).y_mm)); 
        frame_onethird_portion = (0:3).*(movie_length_minute/3)*60*fps; 
        frame_onethird_portion = arrayfun(@(s, e) s+1:min(e, length(tracked_mask)), frame_onethird_portion(1:3), frame_onethird_portion(2:4), 'UniformOutput', false); 
%         frame_10minute = {1:10*60*fps, (10*60*fps+1):(20*60*fps), (20*60*fps+1):min(30*60*fps, length(tracked_mask))};
        mask_onethird_portion = zeros(length(tracked_mask), 3);
        mask_onethird_portion = num2cell(mask_onethird_portion, 1); 
        for i=1:length(mask_onethird_portion)
            mask_onethird_portion{i}(frame_onethird_portion{i}) = 1; 
            mask_onethird_portion{i} = bitand(tracked_mask, mask_onethird_portion{i}); 
        end
        frame_minute_start_end = [1, 60*fps] + repmat(60*fps*(0:(movie_length_minute-1))', 1, 2);
        frame_minute_start_end(end, end) = min(frame_minute_start_end(end, end), length(tracked_mask)); 
        frame_minute = arrayfun(@(s, e) s:e, frame_minute_start_end(:, 1), frame_minute_start_end(:, 2), 'UniformOutput', false); 
        mask_minute = zeros(length(tracked_mask), movie_length_minute);
        mask_minute = num2cell(mask_minute, 1); 
        for i=1:length(mask_minute)
            mask_minute{i}(frame_minute{i}) = 1; 
            mask_minute{i} = bitand(tracked_mask, mask_minute{i}); 
        end
        dist(f) = sum(sqrt(diff(trx(fly_id).x_mm(tracked_mask)).^2 + diff(trx(fly_id).y_mm(tracked_mask)).^2))/1000; 
        dist_onethird_portion(f, :) = cellfun(@(time_mask) sum(sqrt(diff(trx(fly_id).x_mm(logical(time_mask))).^2 + diff(trx(fly_id).y_mm(logical(time_mask))).^2))/1000, ...
            mask_onethird_portion); 
        dist_minute(f, :) = cellfun(@(time_mask) sum(sqrt(diff(trx(fly_id).x_mm(logical(time_mask))).^2 + diff(trx(fly_id).y_mm(logical(time_mask))).^2))/1000, ...
            mask_minute); 
        
        non_orient_mask = labels{f} == 3; 
        if isempty(time_period)
            non_orient_mask = bitand(tracked_mask, reshape(non_orient_mask, size(tracked_mask))); 
            [cluster_labels_cell_by_bout, cluster_bout_length_cell, cluster_labels_changepoint_cell] = convert_cluster_labels_to_bout({non_orient_mask}, 0); 
            dist_in_non_orient(f) = sum(arrayfun(@(s, l) sum(sqrt(diff(trx(fly_id).x_mm(s:s+l-1)).^2 + diff(trx(fly_id).y_mm(s:s+l-1)).^2)), ...
                cluster_labels_changepoint_cell{1}(cluster_labels_cell_by_bout{1}), cluster_bout_length_cell{1}(cluster_labels_cell_by_bout{1})));
            time_in_non_orient(f) = sum(cluster_bout_length_cell{1}(cluster_labels_cell_by_bout{1})/fps); 
        else 
    
            for t=1:length(time_period)
                time_mask = false(size(non_orient_mask));
                time_mask(time_period_start(t):min(time_period_end(t), length(non_orient_mask))) = true;
                non_orient_mask_in_time_period = bitand(non_orient_mask, time_mask); 
                non_orient_mask_in_time_period = bitand(tracked_mask, reshape(non_orient_mask_in_time_period, size(tracked_mask))); 
                [cluster_labels_cell_by_bout, cluster_bout_length_cell, cluster_labels_changepoint_cell] = convert_cluster_labels_to_bout({non_orient_mask_in_time_period}, 0); 
                dist_in_non_orient(f, t) = sum(arrayfun(@(s, l) sum(sqrt(diff(trx(fly_id).x_mm(s:s+l-1)).^2 + diff(trx(fly_id).y_mm(s:s+l-1)).^2)), ...
                    cluster_labels_changepoint_cell{1}(cluster_labels_cell_by_bout{1}), cluster_bout_length_cell{1}(cluster_labels_cell_by_bout{1})));
                time_in_non_orient(f, t) = sum(cluster_bout_length_cell{1}(cluster_labels_cell_by_bout{1})/fps); 
            end
        end
    end
    if same_genotype_pairing
        dist = sum(reshape(dist, 2, []));
        dist_onethird_portion = squeeze(sum(reshape(dist_onethird_portion, 2, [], 3), 1));
        dist_minute = squeeze(sum(reshape(dist_minute, 2, [], 30), 1)); 
        if isempty(time_period)
            speed_in_non_orient = sum(reshape(dist_in_non_orient, 2, []))./sum(reshape(time_in_non_orient, 2, [])); 
        else
            speed_in_non_orient = squeeze(sum(reshape(dist_in_non_orient, 2, [], length(time_period)))./sum(reshape(time_in_non_orient, 2, [], length(time_period)))); 
        end
        if ismember('L_bouts', fieldnames(flymatAll))
            fname = 'L_bouts';
        else
            fname = 'L_total';
        end

        lunge_per_meter = sum(reshape([flymatAll(:).(fname)], 2, []))./dist; 
    else
        if isempty(time_period)
            speed_in_non_orient = dist_in_non_orient./time_in_non_orient; 
        else
            speed_in_non_orient = dist_in_non_orient./time_in_non_orient; 
        end
        lunge_per_meter = [flymatAll(:).(fname)]./dist; 
    end
    
    genotype_list = [flymatAll.genotype]; 
    if same_genotype_pairing
        genotype_list = reshape(genotype_list, 2, []); 
        try 
            assert(all(genotype_list(1, :) == genotype_list(2, :))); 
        catch
            if strcmp(exp_list{e}, 'Tdc2-GAL4_UAS-IR-Irk2_Cyp9f2_CG5808_GM_Aug2019') && find(genotype_list(1, :) ~= genotype_list(2, :)) == 237
                genotype_list(2, 237) = 7; 
            end
        end
        genotype_list = genotype_list(1, :); 
    end
    
    [genotype_mask, genotype_group] = ismember([flymatAll.genotype], genotype_select{e}); 
    dist_cell = arrayfun(@(genotype) dist(genotype_list == genotype), genotype_select{e}, 'UniformOutput', false); 
    lunge_per_meter_per_pair_cell = arrayfun(@(genotype) lunge_per_meter(genotype_list == genotype), genotype_select{e}, 'UniformOutput', false); 
    dist_onethird_portion_cell = arrayfun(@(genotype) dist_onethird_portion(genotype_list == genotype, :), genotype_select{e}, 'UniformOutput', false); 
    dist_minute_cell = arrayfun(@(genotype) dist_minute(genotype_list == genotype, :), genotype_select{e}, 'UniformOutput', false); 
    speed_in_non_orient_cell = arrayfun(@(genotype) speed_in_non_orient(genotype_list == genotype), genotype_select{e}, 'UniformOutput', false); 
    dist_struct = cell2struct(dist_cell, strcat('genotype', strip(cellstr(num2str([genotype_select{e}]')))), 2); 
    lunge_per_meter_per_pair_struct = cell2struct(lunge_per_meter_per_pair_cell, strcat('genotype', strip(cellstr(num2str([genotype_select{e}]')))), 2); 
    dist_onethird_portion_struct = cell2struct(dist_onethird_portion_cell, strcat('genotype', strip(cellstr(num2str([genotype_select{e}]')))), 2); 
    dist_minute_struct = cell2struct(dist_minute_cell, strcat('genotype', strip(cellstr(num2str([genotype_select{e}]')))), 2); 
    speed_in_non_orient_struct = cell2struct(speed_in_non_orient_cell, strcat('genotype', strip(cellstr(num2str([genotype_select{e}]')))), 2); 
    dist_travel_struct(e) = struct('exp', exp_list{e}, sprintf('dist_%s', suffix), dist_struct, sprintf('dist_%dminute_%s', round(movie_length_minute/3), suffix), dist_onethird_portion_struct, sprintf('dist_minute_%s', suffix), dist_minute_struct, sprintf('lungePerMeter_%s', suffix), lunge_per_meter_per_pair_struct, sprintf('speed_in_non_orient_%s', suffix), speed_in_non_orient_struct); 
end

cd(output_path); 
save('dist_travel-lunge_per_meter.mat', 'dist_travel_struct'); 

%% box plot
% a fairly special case for plotting speed_in_non_orient in separate time
% periods for pre-opto, opto and post-opto
figure()
genotype_select = repmat(genotype_select, length(time_period), 1);
genotype_select = cell2mat(genotype_select) + 0.1*[1:length(time_period)]'; 
genotype_select = num2cell(genotype_select(:)); 

genotype_name_time_period = cellfun(@(g) cellfun(@(t) sprintf('%s-%s', g, t), time_period, 'UniformOutput', false), genotype_name{1}, 'UniformOutput', false);
genotype_name_time_period = horzcat(genotype_name_time_period{:}); 

speed_in_non_orient_temp = cellfun(@(d) mat2cell(d, size(d, 1), ones(1, length(time_period))), speed_in_non_orient_cell, 'UniformOutput', false); 
speed_in_non_orient_temp = horzcat(speed_in_non_orient_temp{:}); 
boxplot(vertcat(speed_in_non_orient_temp{:}), repelem(1:length(genotype_select), cellfun(@length, speed_in_non_orient_temp)'), ...
    'Symbol', '', 'Labels', genotype_name_time_period);
hold on;
scatter(repelem(1:length(genotype_select), cellfun(@length, speed_in_non_orient_temp)'), vertcat(speed_in_non_orient_temp{:}), 12); 
ylabel('speed during non-orient'); 
xtickangle(45);
hold off; 
set(gcf, 'Renderer', 'painters'); 
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, fullfile(output_path, sprintf('%s-speed_during_non_orient.eps', exp_list{1})), 'epsc');
saveas(gcf, fullfile(output_path, sprintf('%s-speed_during_non_orient.png', exp_list{1})));

% for e=1:length(exp_list)
%     exp_path = fullfile(base_path, exp_list{e}); 
%     cd(exp_path);
%     flymat_file = dir('FLYMAT_*.mat');
%     load(flymat_file.name); %strcat('FLYMAT_', exp_list{e})); 
%    
%     [genotype_mask, genotype_group] = ismember([flymatAll.genotype], genotype_select{e}); 
%     dist = dist_travel_struct(e).dist(genotype_mask); 
%     genotype = [flymatAll(genotype_mask).genotype]; 
%     [p_val_kruskal_wallis, ~, stats] = kruskalwallis(dist, genotype, 'off'); 
%     c = multcompare(stats, 'CType', 'bonferroni', 'Display', 'off'); 
%     figure();
%     boxplot(dist, genotype, ...
%         'Symbol', '', 'Labels', cellfun(@(name, num) sprintf('%s; genotype %d', name, num), genotype_name{e}, num2cell(genotype_select{e}), 'UniformOutput', false));
%     hold on;
%     scatter(genotype_group(genotype_group>0), dist, 12); 
%     ylabel('total distance traveled (m)'); 
%     stat_report_str = 'multi-comparison (bonferroni-corrected): ';
%     for p=1:size(c, 1)
%         stat_report_str = strcat(stat_report_str, sprintf('%s-v-%s: %.5f; ', stats.gnames{c(p, 1)}, stats.gnames{c(p, 2)}, c(p, end))); 
%     end
%     title({sprintf('Kruskal-Wallis %.5f', p_val_kruskal_wallis), stat_report_str}); 
%     xtickangle(45);
%     hold off; 
%     set(gcf, 'Renderer', 'painters'); 
%     set(gcf, 'Position', get(0, 'Screensize'));
%     saveas(gcf, fullfile(output_path, sprintf('%s-distance_traveled.eps', exp_list{e})), 'epsc');
%     saveas(gcf, fullfile(output_path, sprintf('%s-distance_traveled.png', exp_list{e})));
    
%     speed = dist_travel_struct(e).speed(genotype_mask); 
%     [p_val_kruskal_wallis, ~, stats] = kruskalwallis(speed, genotype, 'off'); 
%     c = multcompare(stats, 'CType', 'bonferroni', 'Display', 'off'); 
%     figure();
%     boxplot(speed, genotype, ...
%         'Symbol', '', 'Labels', cellfun(@(name, num) sprintf('%s; genotype %d', name, num), genotype_name{e}, num2cell(genotype_select{e}), 'UniformOutput', false));
%     hold on;
%     scatter(genotype_group(genotype_group>0), speed, 12); 
%     ylabel('mean speed in non-orient bouts (mm/s)'); 
%     stat_report_str = 'multi-comparison (bonferroni-corrected): ';
%     for p=1:size(c, 1)
%         stat_report_str = strcat(stat_report_str, sprintf('%s-v-%s: %.5f; ', stats.gnames{c(p, 1)}, stats.gnames{c(p, 2)}, c(p, end))); 
%     end
%     title({sprintf('Kruskal-Wallis %.5f', p_val_kruskal_wallis), stat_report_str}); 
%     xtickangle(45);
%     hold off; 
%     set(gcf, 'Renderer', 'painters'); 
%     set(gcf, 'Position', get(0, 'Screensize'));
%     saveas(gcf, fullfile(output_path, sprintf('%s-mean_speed.eps', exp_list{e})), 'epsc');
%     saveas(gcf, fullfile(output_path, sprintf('%s-mean_speed.png', exp_list{e})));
% end