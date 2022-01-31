output_folder = 'D:\xubo\nervy_w_Kenichi\figures-122921\box_plot-genetic_manipulations'; 

% ===== WT and KO =====
exp_name = 'WT_and_KO'; 
base_path = 'D:\xubo\nervy_w_Kenichi\figures-121321\behav_bout_count'; 
count_file = {'frame_count-CG3385KO_CS-bc-11_Df-Exel6082_CS-bc-6_GM_Aug2019.mat', 'frame_count-CG3385_KO_CS-bc-11_SM_Mar2017.mat'};
genotype_select = {[4, 1], [3, 99]};
genotype_name = {{'delta-nvy GH', 'wildtype GH'}, {'delta-nvy SH', 'wildtype SH'}};
stat_compare = true; 
pairwise_comp_rel_idx = [1, 6, 2, 5]; 

% ==== RNAi =====
% exp_name = 'rnai';
% base_path = 'D:\xubo\nervy_w_Kenichi\figures-122921\box_plot-genetic_manipulations';
% count_file = {'frame_count-Elav-GAL4_UAS-IR-nvy(VDRC100373,JF03349)_May2018.mat'};
% genotype_select = {[1, 2, 3]};
% genotype_name = {{'Elav-GAL4, 260b', '+, UAS-IR nvy', 'Elav-GAL4, UAS-IR nvy'}}; 
% stat_compare = false; 

% ===== rescue GH =====
% exp_name = 'rescue'; 
% base_path = 'D:\xubo\nervy_w_Kenichi\figures-122921\box_plot-genetic_manipulations';
% count_file = {'frame_count-CG3385KO(CS-bc-11)_Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_GM_Jul2017.mat'};
% genotype_select = {[2, 3, 1]};  
% genotype_name = {{'delta-nvy, elav-GAL4>GFP', 'delta-nvy, +>nvy', 'delta-nvy, elav-GAL4>nvy'}};
% stat_compare = false; 

% ===== over-expression SH =====
% exp_name = 'over_expression'; 
% base_path = 'D:\xubo\nervy_w_Kenichi\figures-122921\box_plot-genetic_manipulations';
% count_file = {'frame_count-Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_SM_Aug2017'};
% genotype_select = {[21, 31, 11]}; 
% genotype_name = {{'Elav-GAL4, UAS-GFP', '+, UAS-nvy', 'Elav-GAL4, UAS-nvy'}};
% stat_compare = false; 

% ===== nvy KO paired with mated female =====
% exp_name = 'male_v_female';
% base_path = 'D:\xubo\nervy_w_Kenichi\figures-122921\box_plot-genetic_manipulations';
% count_file = {'behav_count-CG3385_KO_CS-bc-11_GMvsWT_GMF_60min_May2017-2.mat'}; 
% genotype_select = {[9, 1]}; 
% genotype_name = {{'wildtype GH male', 'delta-nvy GH male'}};
% stat_compare = false; 

count_field = 'frame_count_pair'; 
for e=1:length(count_file)
    temp = load(fullfile(base_path, count_file{e}), count_field); 
    if e == 1
        count = struct2cell(temp.(count_field))'; 
        fieldname = fieldnames(temp.(count_field)); 
    else
        count = [count, struct2cell(temp.(count_field))']; 
        fieldname = [fieldname; fieldnames(temp.(count_field))]; 
    end
end
genotype_select = horzcat(genotype_select{:}); 
genotype_name = horzcat(genotype_name{:}); 

[~, idx] = ismember(fieldname, arrayfun(@(n) sprintf('genotype%d', n), genotype_select, 'UniformOutput', false)); 
count = count(:, idx); 

num_of_label = 6; 
behav_name = {'stopping', 'orienting', 'non-orienting', 'lunging', 'wing extension', 'fragment'};

for b=1:num_of_label
    count_vec = horzcat(count{b, :}); 
    category_vec = repelem(1:length(genotype_select), cellfun(@length, count(b, :))); 
    
    if stat_compare
        [p_val_kruskal_wallis, ~, stats] = kruskalwallis(count_vec, category_vec, 'off'); 
        c = multcompare(stats, 'CType', 'bonferroni', 'Display', 'on'); 
        p_val_multicomp = c(pairwise_comp_rel_idx, end); 
%         for p=1:length(pairwise_comp_rel_idx)
%             [pval, ~, rks_stats] = ranksum([count{b, c(pairwise_comp_rel_idx(p), 1)}], [count{b, c(pairwise_comp_rel_idx(p), 2)}]); 
%             fprintf("%d %d %.6f %f\n", c(pairwise_comp_rel_idx(p), 1), c(pairwise_comp_rel_idx(p), 2), pval, rks_stats.ranksum - stats.n(c(pairwise_comp_rel_idx(p), 1))); 
%         end
    end
    pct = cellfun(@(c) prctile(c, [25, 50, 75]), count(b, :), 'UniformOutput', false); 
    pct_all = vertcat(pct{:});
    
    figure();
    boxplot(horzcat(count{b, :}), category_vec', 'Symbol', '', 'Labels', genotype_name);
    hold on;
    scatter(category_vec', count_vec, 12); 
    field_part = strsplit(count_field, '_'); 
    ylabel(sprintf('# %ss %s/%s', field_part{1}, behav_name{b}, field_part{3})); 
    xtickangle(45);
    hold off; 
    if stat_compare
        title({sprintf('Kruskal-Wallis %.5f', p_val_kruskal_wallis),...
            sprintf('multi-comparison (bonferroni-corrected) \n dnvyGH-v-WTGH: %.5f; dnvySH-v-WTSH: %.5f; \n dnvyGH-v-dnvySH: %.5f; WTGH-v-WTSH:  %.5f', ...
            p_val_multicomp)}); 
    end
    set(gcf, 'Renderer', 'painters'); 
    set(gcf, 'Position', get(0, 'Screensize'));
%     saveas(gcf, fullfile(output_folder, sprintf('%s-%s_%s.eps', exp_name, behav_name{b}, count_field)), 'epsc');
%     saveas(gcf, fullfile(output_folder, sprintf('%s-%s_%s.png', exp_name, behav_name{b}, count_field)));
end