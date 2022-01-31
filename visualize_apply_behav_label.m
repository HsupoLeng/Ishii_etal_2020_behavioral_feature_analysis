output_folder = 'D:\xubo\nervy_w_Kenichi\data-013122'; 

num_of_label = 6;
behav_name = {'stopping', 'orienting', 'non-orienting', 'lunging', 'wing extension', 'fragment'};
remove_fragment = true; 
n_repeat = 10000;
n_tail = 5; 
fps = 60; 
same_genotype_pairing = true; 

% ===== KO =====
base_path = 'D:\xubo\behavioral_movies\CG3385KO_CS-bc-11_Df-Exel6082_CS-bc-6_GM_Aug2019'; 
flymat_file = 'FLYMAT_CG3385KO_CS-bc-11_Df-Exel6082_CS-bc-6_GM_Aug2019.mat'; 
genotype_select = {1, 4}; 
genotype_name = {'wildtype GH', 'delta-nvy GH'};
time_period = {}; 

% base_path = {'D:\xubo\behavioral_movies\CG3385KO_CS-bc-11_Df-Exel6082_CS-bc-6_GM_Aug2019', 'Z:\Kenichi\CG3385_KO_CS-bc-11_SM_Mar2017'};
% flymat_file = {'FLYMAT_CG3385KO_CS-bc-11_Df-Exel6082_CS-bc-6_GM_Aug2019.mat', 'FLYMAT_CG3385_KO_CS-bc-11_SM_Mar2017.mat'};
% genotype_select = {1, 99};
% genotype_name = {'wildtype GH', 'wildtype SH'};
% time_period = {}; 

% base_path = {'D:\xubo\behavioral_movies\CG3385KO_CS-bc-11_Df-Exel6082_CS-bc-6_GM_Aug2019', 'Z:\Kenichi\CG3385_KO_CS-bc-11_SM_Mar2017'};
% flymat_file = {'FLYMAT_CG3385KO_CS-bc-11_Df-Exel6082_CS-bc-6_GM_Aug2019.mat', 'FLYMAT_CG3385_KO_CS-bc-11_SM_Mar2017.mat'};
% genotype_select = {4, 3};
% genotype_name = {'delta-nvy GH', 'delta-nvy SH'};
% time_period = {}; 

% base_path = 'Z:\Kenichi\CG3385_KO_CS-bc-11_SM_Mar2017';
% flymat_file = 'FLYMAT_CG3385_KO_CS-bc-11_SM_Mar2017.mat';
% genotype_select = {99, 3}; 
% genotype_name = {'wildtype SH', 'delta-nvy SH'};
% time_period = {}; 

% ==== RNAi =====
% base_path = 'Z:\Kenichi\Elav-GAL4_UAS-IR-nvy(VDRC100373,JF03349)_May2018'; 
% flymat_file = 'FLYMAT_Elav-GAL4_UAS-IR-nvy(VDRC100373,JF03349)_May2018.mat'; 
% genotype_select = {1, 2};
% genotype_name = {'+(260b) / + ; + / Elav-GAL4', 'UAS-IR-nvy(VDRC100273) / +(CS) ; +'}; 
% time_period = {}; 

% base_path = 'Z:\Kenichi\Elav-GAL4_UAS-IR-nvy(VDRC100373,JF03349)_May2018'; 
% flymat_file = 'FLYMAT_Elav-GAL4_UAS-IR-nvy(VDRC100373,JF03349)_May2018.mat'; 
% genotype_select = {1, 3};
% genotype_name = {'Elav-GAL4, 260b', 'Elav-GAL4, UAS-IR nvy'}; 
% time_period = {}; 

% base_path = 'Z:\Kenichi\Elav-GAL4_UAS-IR-nvy(VDRC100373,JF03349)_May2018'; 
% flymat_file = 'FLYMAT_Elav-GAL4_UAS-IR-nvy(VDRC100373,JF03349)_May2018.mat'; 
% genotype_select = {2, 3};
% genotype_name = {'+, UAS-IR nvy', 'Elav-GAL4, UAS-IR nvy'}; 
% time_period = {}; 

% base_path = 'Z:\Kenichi\Elav-GAL4_UAS-IR-nvy(VDRC100373,JF03349)_May2018'; 
% flymat_file = 'FLYMAT_Elav-GAL4_UAS-IR-nvy(VDRC100373,JF03349)_May2018.mat'; 
% genotype_select = {[1, 2], 3};
% genotype_name = {'controls', 'UAS-IR-nvy(VDRC100273) / + ; + / Elav-GAL4'}; 
% genotype_select = {1, 2, 3}; 
% genotype_name = {'Elav-GAL4, 260b', '+, UAS-IR nvy', 'Elav-GAL4, UAS-IR nvy'}; 
% time_period = {}; 

% ===== over-expression SH =====
% base_path = 'Z:\Kenichi\Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_SM_Aug2017'; 
% flymat_file = 'FLYMAT_Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_SM_Aug2017'; 
% genotype_select = {21, 31, 11};
% genotype_name = {'elav-GAL4>GFP', '+>nvy', 'elav-GAL4>nvy'}; 
% time_period = {}; 

% base_path = 'Z:\Kenichi\Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_SM_Aug2017'; 
% flymat_file = 'FLYMAT_Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_SM_Aug2017'; 
% genotype_select = {31, 11};
% genotype_name = {'+, UAS-nvy', 'Elav-GAL4, UAS-nvy'}; 
% time_period = {}; 
 
% base_path = 'Z:\Kenichi\Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_SM_Aug2017'; 
% flymat_file = 'FLYMAT_Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_SM_Aug2017'; 
% genotype_select = {21, 11};
% genotype_name = {'Elav-GAL4, UAS-GFP', 'Elav-GAL4, UAS-nvy)'}; 
% time_period = {}; 

% base_path = 'Z:\Kenichi\Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_SM_Aug2017'; 
% flymat_file = 'FLYMAT_Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_SM_Aug2017'; 
% genotype_select = {21, 31};
% genotype_name = {'elav-GAL4 / 10xUAS-IVS-Syn21-GFP-p10 (VK00005, CS-bc-6)', '+ / 10xUAS-IVS-myc-nvy (VK00005, CS-bc-6)'}; 
% time_period = {}; 

% ===== rescue GH =====
% base_path = 'Z:\Kenichi\CG3385KO(CS-bc-11)_Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_GM_Jul2017';
% flymat_file = 'FLYMAT_CG3385KO(CS-bc-11)_Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_GM_Jul2017.mat';
% genotype_select = {2, 3, 1}; 
% genotype_name = {'delta-nvy, elav-GAL4>GFP', 'delta-nvy, +>nvy', 'delta-nvy, elav-GAL4>nvy'};
% time_period = {}; 

% base_path = 'Z:\Kenichi\CG3385KO(CS-bc-11)_Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_GM_Jul2017';
% flymat_file = 'FLYMAT_CG3385KO(CS-bc-11)_Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_GM_Jul2017.mat';
% genotype_select = {2, 1}; 
% genotype_name = {'delta-nvy, Elav-GAL4, UAS-GFP', 'delta-nvy, elav-gal4, UAS-nvy'};
% time_period = {}; 

% base_path = 'Z:\Kenichi\CG3385KO(CS-bc-11)_Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_GM_Jul2017';
% flymat_file = 'FLYMAT_CG3385KO(CS-bc-11)_Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_GM_Jul2017.mat';
% genotype_select = {3, 1}; 
% genotype_name = {'delta-nvy, +, UAS-nvy', 'delta-nvy, Elav-GAL4, UAS-nvy'};
% time_period = {}; 
 
% base_path = 'Z:\Kenichi\CG3385KO(CS-bc-11)_Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_GM_Jul2017';
% flymat_file = 'FLYMAT_CG3385KO(CS-bc-11)_Elav-orNo-gal4,UAS-nvyorGFP(CS-bc-6)_GM_Jul2017.mat';
% genotype_select = {2, 3}; 
% genotype_name = {'delta-nvy, elav>GFP', 'delta-nvy, +>nvy'};
% time_period = {}; 

% ===== Heisenberg =====
% base_path = 'Z:\Kenichi\CSMH_SM_vsWT_GM_HeisenbergChamber_May2019';
% flymat_file = 'FLYMAT_CSMH_SM_vsWT_GM_HeisenbergChamber_May2019.mat';
% genotype_select = {90, 91}; 
% genotype_name = {'wildtype GH', 'wildtype SH'}; 
% time_period = {}; 

% base_path = 'Z:\Kenichi\CG3385_KO_CS-bc-11_GM_vsWT_GM_HeisenbergChamber_May2019';
% flymat_file = 'FLYMAT_CG3385_KO_CS-bc-11_GM_vsWT_GM_HeisenbergChamber_May2019.mat';
% genotype_select = {90, 1};
% genotype_name = {'wildtype GH', 'delta-nvy GH'};
% time_period = {}; 

% base_path = 'Z:\Kenichi\CG3385_KO_CS-bc-11_SMvsCSMH_GM_HeisenbergChamber_Aug2021';
% flymat_file = 'FLYMAT_CG3385_KO_CS-bc-11_SMvsCSMH_GM_HeisenbergChamber_Aug2021.mat';
% genotype_select = {99, 11};
% genotype_name = {'wildtype GH', 'delta-nvy SH'}; 
% time_period = {}; 

% ===== nvy KO paired with mated female =====
% base_path = 'Z:\Kenichi\CG3385_KO_CS-bc-11_GMvsWT_GMF_60min_May2017-2'; 
% flymat_file = 'FLYMAT_CG3385_KO_CS-bc-11_GMvsWT_GMF_60min_May2017-2_newHBclassifier.mat'; 
% genotype_select = {9, 1};
% genotype_name = {'wildtype GH male', 'delta-nvy GH male'};
% time_period = {}; 

% ===== nvy+/Tdc2+ neurons silencing GM =====
% base_path = 'Z:\Kenichi\Tdc2-GAL4_nvy-LexA_Tub-FRT-GAL80-FRTorTub-FRT-stop-FRT-GAL80_LexAop-FLPL_UAS-Kir_GM_Nov2018'; 
% flymat_file = 'FLYMAT_Tdc2-GAL4,nvy-LexA_Tub-FRT-GAL80-FRTorTub-FRT-stop-FRT-GAL80,LexAop-FLPL_UAS-Kir_GM_Nov2018.mat'; 
% genotype_select = {23, 43, 63};
% genotype_name = {'Tdc2-GAL4, nvy-LexA, Tub>GAL80>, UAS-Kir2.1', 'Tdc2-GAL4, nvy-LexA, Tub>GAL80>, LexA-FLPL, UAS-GFP', 'Tdc2-GAL4, nvy-LexA, Tub>GAL80>, LexA-FLPL, UAS-Kir2.1'}; 
% time_period = {}; 

% ===== nvy+/Tdc2+ neurons silencing SM =====
% base_path = 'Z:\Kenichi\Tdc2-GAL4_nvy-LexA_Tub-FRT-GAL80-FRTorTub-FRT-stop-FRT-GAL80_LexAop-FLPL_UAS-Kir_SM_Nov2018';
% flymat_file = 'FLYMAT_Tdc2-GAL4,nvy-LexA_Tub-FRT-GAL80-FRTorTub-FRT-stop-FRT-GAL80,LexAop-FLPL_UAS-Kir_SM_Nov2018.mat';
% genotype_select = {11, 21, 41, 61};
% genotype_name = {'Tdc2-GAL4, nvy-LexA, UAS-Kir2.1', 'Tdc2-GAL4, nvy-LexA, Tub>GAL80>, UAS-Kir2.1', 'Tdc2-GAL4, nvy-LexA, Tub>GAL80>, LexA-FLPL, UAS-GFP', 'Tdc2-GAL4, nvy-LexA, Tub>GAL80>, LexA-FLPL, UAS-Kir2.1'}; 
% time_period = {}; 

% ===== nvy-/Tdc2+ neurons silencing GM =====
% base_path = 'Z:\Kenichi\Tdc2-GAL4_nvy-LexA_Tub-FRT-GAL80-FRTorTub-FRT-stop-FRT-GAL80_LexAop-FLPL_UAS-Kir_GM_Nov2018'; 
% flymat_file = 'FLYMAT_Tdc2-GAL4,nvy-LexA_Tub-FRT-GAL80-FRTorTub-FRT-stop-FRT-GAL80,LexAop-FLPL_UAS-Kir_GM_Nov2018.mat'; 
% genotype_select = {33, 53, 73};
% genotype_name = {'Tdc2-GAL4, nvy-LexA, Tub>stop>GAL80, UAS-Kir2.1', 'Tdc2-GAL4, nvy-LexA, Tub>stop>GAL80, LexA-FLPL, UAS-GFP', 'Tdc2-GAL4, nvy-LexA, Tub>stop>GAL80, LexA-FLPL, UAS-Kir2.1'}; 
% time_period = {}; 

% ===== nvy+/Tdc2+ neurons silencing SM =====
% base_path = 'Z:\Kenichi\Tdc2-GAL4_nvy-LexA_Tub-FRT-GAL80-FRTorTub-FRT-stop-FRT-GAL80_LexAop-FLPL_UAS-Kir_SM_Nov2018';
% flymat_file = 'FLYMAT_Tdc2-GAL4,nvy-LexA_Tub-FRT-GAL80-FRTorTub-FRT-stop-FRT-GAL80,LexAop-FLPL_UAS-Kir_SM_Nov2018.mat';
% genotype_select = {31, 51, 71};
% genotype_name = {'Tdc2-GAL4, nvy-LexA, Tub>stop>GAL80, UAS-Kir2.1', 'Tdc2-GAL4, nvy-LexA, Tub>stop>GAL80, LexA-FLPL, UAS-GFP', 'Tdc2-GAL4, nvy-LexA, Tub>stop>GAL80, LexA-FLPL, UAS-Kir2.1'}; 
% time_period = {}; 

% ===== Optogenetic activation of nvy+/Tdc2+ =====
% base_path = 'Z:\Kenichi\Tdc2-GAL4_nvy-LexA_Tub-FRT-GAL80_LexAop-FLPL_UAS-CsChrim_SM_vsCSmale_LED5minOffOnOff_Feb2019'; 
% flymat_file = 'FLYMAT_Tdc2-GAL4_nvy-LexA_Tub-FRT-GAL80_LexAop-FLPL_UAS-CsChrim_SM_vsCSmale_LED5minOffOnOff_Feb2019'; 
% genotype_select = {11, 21, 31};
% genotype_name = {'nvy-LexA, Tdc2-GAL4, Tub>GAL80>, +, UAS-CsChrimson', 'nvy-LexA, Tdc2-GAL4, Tub>GAL80>, LexAop2-FLPL, UAS-GFP', 'nvy-LexA, Tdc2-GAL4, Tub>GAL80>, LexAOp2-FLPL, UAS-CsChrimson'}; 
% time_period = {'pre-opto', 'opto', 'post-opto'};
% time_period_length = 5*60*fps; 

if isempty(time_period)
    [info, cluster_labels_cell_by_bout, cluster_bout_length_cell, cluster_labels_changepoint_cell, labels_all] = ...
        apply_behav_label(base_path, flymat_file, genotype_select, genotype_name, remove_fragment);
else
    cluster_labels_cell_by_bout = cell(length(time_period), length(genotype_select)); 
    cluster_bout_length_cell = cell(size(cluster_labels_cell_by_bout)); 
    cluster_labels_changepoint_cell = cell(size(cluster_labels_cell_by_bout));
    
    time_period_start = [0:length(time_period)-1]' .* time_period_length + 1;
    time_period_end = time_period_start + time_period_length - 1; 
     
    remove_fragment_temp = false;
    [~, ~, ~, ~, labels_all] = ...
        apply_behav_label(base_path, flymat_file, genotype_select, genotype_name, remove_fragment_temp);
    for g=1:length(genotype_select)
        for t=1:length(time_period)
            labels_in_period = cellfun(@(f) f(time_period_start(t):min(time_period_end(t), length(f))), labels_all{g}, 'UniformOutput', false); 
      
            if remove_fragment
                for f=1:length(labels_in_period)
                    labels_one_fly = labels_in_period{f}; 
                    labels_in_period{f} = labels_one_fly(labels_one_fly ~= num_of_label); 
                end
            end
        
            ind = sub2ind([length(time_period), length(genotype_select)], t, g); 
            [cluster_labels_cell_by_bout{ind}, cluster_bout_length_cell{ind}, cluster_labels_changepoint_cell{ind}] = convert_cluster_labels_to_bout(labels_in_period, 0); 
        end
    end
    genotype_select = repmat(genotype_select, length(time_period), 1);
    genotype_select = cell2mat(genotype_select) + 0.1*[1:length(time_period)]'; 
    genotype_select = num2cell(genotype_select(:)); 
    
    genotype_name_time_period = cellfun(@(g) cellfun(@(t) sprintf('%s-%s', g, t), time_period, 'UniformOutput', false), genotype_name, 'UniformOutput', false); 
    genotype_name_time_period = horzcat(genotype_name_time_period{:}); 
end
%% visualize results
% Visuailze bout length histogram for each behavior
% ax = [];
% for b=1:num_of_label
%     figure();
%     if ismember(b, [1, 4])
%         bin_w = 3; 
%     else
%         bin_w = 3;
%     end
%     for g=1:length(genotype_select)
%         ax(g) = subplot(1, 2, g);
%         behav_mask = cellfun(@(l) l == b, cluster_labels_cell_by_bout{g}, 'UniformOutput', false);
%         bout_length_one_behav = cellfun(@(l, ma) l(ma), cluster_bout_length_cell{g}, behav_mask, 'UniformOutput', false); 
%         histogram(vertcat(bout_length_one_behav{:}), 'BinWidth', 1);
%         title(genotype_name{g}); 
%         set(gca, 'YScale', 'log');
%         xlabel('bout length (# of frames)'); 
% 
%         if g == 1
%             ylabel('# of bouts'); 
%         end
% 
%         num_short_bout = cellfun(@(l) sum(l <= 3), bout_length_one_behav, 'UniformOutput', false);
%         ratio_short_bout = sum([num_short_bout{:}])/sum(cellfun(@length, bout_length_one_behav));
%         fprintf('%.3f%% bouts of behavior %s for %s is shorter than threshold\n', ratio_short_bout*100, behav_name{b}, genotype_name{g});  
%     end
%     suptitle(behav_name{b});
%     linkaxes(ax, 'y');
% end
% 

% % Visualize bout interval histogram for each behavior
% bout_intv_cell = cell(1, 2); %cell(num_of_label, 2);
% for g=1:length(genotype_select)
%     for b=4 %1:num_of_label
%         bout_idx = cellfun(@(l) find(l == b), cluster_labels_cell_by_bout{g}, 'UniformOutput', false); 
%         bout_mask_one_apart = cellfun(@(b_idx) diff(b_idx) == 2, bout_idx, 'UniformOutput', false); 
%         bout_idx_one_apart = cellfun(@(b_idx, mask) b_idx([false; mask]), bout_idx, bout_mask_one_apart, 'UniformOutput', false); 
%         bout_inbetween = cellfun(@(b_idx) b_idx - 1, bout_idx_one_apart, 'UniformOutput', false);
%         bout_intv_cell{g} = cellfun(@(b_idx, b_length) b_length(b_idx), bout_inbetween, cluster_bout_length_cell{g}, 'UniformOutput', false);
%     end
% end
% ax = [];
% % for b=1:num_of_label
% figure();
% bin_w = 3; 
% for g=1:length(genotype_select)
%     ax(g) = subplot(1, 2, g);
%     bout_intv_length_one_behav = bout_intv_cell{g}; 
%     histogram(vertcat(bout_intv_length_one_behav{:}), 'BinWidth', 1, 'BinLimits', [0, 50]);
%     title(genotype_name{g}); 
%     set(gca, 'YScale', 'log');
%     xlabel('bout interval length (# of frames)'); 
% 
%     if g == 1
%         ylabel('# of bouts'); 
%     end
% end
% suptitle({'Interval between lunge bouts'});
% linkaxes(ax, 'y');
% end

%% visualize bout label and length histograms over all behaviors
% cluster_labels_cell_rec = cell(size(cluster_labels_cell_by_bout)); 
% for g=1:length(genotype_select)
%     cluster_labels_cell_rec{g} = cellfun(@(bo, le) repelem(bo, le), cluster_labels_cell_by_bout{g}, cluster_bout_length_cell{g}, 'UniformOutput', false); 
% end
% figure(); 
% for g=1:length(genotype_select)
%     ax(g) = subplot(1, 2, g);
%     histogram(vertcat(cluster_labels_cell_rec{g}{:}), 'BinLimits', [0.5, 6.5], 'BinWidth', 1); 
%     set(gca, 'YScale', 'log'); 
%     title(genotype_name{g});
% end
% linkaxes(ax, 'y');
% figure(); 
% for g=1:length(genotype_select)
%     ax(g) = subplot(1, 2, g);
%     histogram(vertcat(cluster_labels_cell_rec{g}{:}), 'BinLimits', [0.5, 6.5], 'BinWidth', 1, 'Normalization', 'probability');  
%     title(genotype_name{g});
% end
% linkaxes(ax, 'y');
% 
% figure();
% ax = []; 
% for g=1:length(genotype_select)
%     ax(g) = subplot(1, 2, g); 
%     histogram(vertcat(cluster_bout_length_cell{g}{:}), 'BinWidth', 3); 
%     set(gca, 'YScale', 'log'); 
%     title(genotype_name{g});
% end
% linkaxes(ax, 'y');
%  
% figure(); 
% for g=1:length(genotype_select)
%     ax(g) = subplot(1, 2, g);
%     histogram(vertcat(cluster_labels_cell_by_bout{g}{:}), 'BinLimits', [0.5, 6.5], 'BinWidth', 1); 
%     set(gca, 'YScale', 'log');
%     title(genotype_name{g});
% end
% linkaxes(ax, 'y');
% 
% figure(); 
% for g=1:length(genotype_select)
%     ax(g) = subplot(1, 2, g);
%     histogram(vertcat(cluster_labels_cell_by_bout{g}{:}), 'BinLimits', [0.5, 6.5], 'BinWidth', 1, 'Normalization', 'probability'); 
%     title(genotype_name{g});
% end
% linkaxes(ax, 'y');

%% visualize bout count and bout length comparison
if ischar(base_path)
    [~, dataset_name, ~] = fileparts(base_path); 
else
    dataset_name_cell = {};
    for g=1:length(base_path)
        [~, dataset_name_cell{g}, ~] = fileparts(base_path{g}); 
    end
    dataset_name = strjoin(dataset_name_cell, '_vs_'); 
end
genotype_fieldname = strcat('genotype', strip(cellfun(@(s) strrep(s, '.', '_'), cellstr(num2str([genotype_select{:}]')), 'UniformOutput', false))); 

labels_all = cell2struct(labels_all, genotype_fieldname, 2); 
% save(fullfile(output_folder, sprintf('behavior_labels-%s.mat', dataset_name)), 'labels_all', 'genotype_name', 'genotype_select');

%     
% ax = [];
bout_count_fly = cell(num_of_label, length(genotype_select)); 
frame_count_fly = cell(num_of_label, length(genotype_select)); 
bout_count_fly_normalized = cell(size(bout_count_fly));
frame_count_fly_normalized = cell(size(frame_count_fly));

for g=1:length(genotype_select)
    for b=1:num_of_label
         behav_mask = cellfun(@(l) l == b, cluster_labels_cell_by_bout{g}, 'UniformOutput', false);
         bout_count_fly{b, g} = cellfun(@sum, behav_mask); 
         n_bouts_per_fly = cellfun(@length, cluster_labels_cell_by_bout{g}); 
         bout_count_fly_normalized{b, g} = bout_count_fly{b, g}./n_bouts_per_fly; 
        
         frame_count_fly{b, g} = cellfun(@(l, m) sum(l(m)), cluster_bout_length_cell{g}, behav_mask);
         n_frames_per_fly = cellfun(@sum, cluster_bout_length_cell{g}); 
         frame_count_fly_normalized{b, g} = frame_count_fly{b, g}./n_frames_per_fly; 
     end
end

for b=1:num_of_label
    figure()
    if same_genotype_pairing
        boxplot(horzcat(bout_count_fly{b, :}), repelem(1:length(genotype_select), cellfun(@length, bout_count_fly(b, :)))', ...
            'Symbol', '', 'Labels', genotype_name);
    else
        boxplot(horzcat(bout_count_fly{b, :}), repelem(1:length(genotype_select), cellfun(@length, bout_count_fly(b, :)))', ...
            'Symbol', '', 'Labels', genotype_name_time_period);
    end
    hold on;
    scatter(repelem(1:length(genotype_select), cellfun(@length, bout_count_fly(b, :)))', horzcat(bout_count_fly{b, :}), 12); 
    ylabel(sprintf('# bouts %s/fly', behav_name{b})); 
    xtickangle(45);
    hold off; 
    set(gcf, 'Renderer', 'painters'); 
    set(gcf, 'Position', get(0, 'Screensize'));
%     saveas(gcf, fullfile(output_folder, sprintf('%s-%s_bout_count.eps', dataset_name, behav_name{b})), 'epsc');
%     saveas(gcf, fullfile(output_folder, sprintf('%s-%s_bout_count.png', dataset_name, behav_name{b})));
     
    figure()
    if same_genotype_pairing
        boxplot(horzcat(frame_count_fly{b, :}), repelem(1:length(genotype_select), cellfun(@length, frame_count_fly(b, :)))', ...
            'Symbol', '', 'Labels', genotype_name);
    else
        boxplot(horzcat(frame_count_fly{b, :}), repelem(1:length(genotype_select), cellfun(@length, frame_count_fly(b, :)))', ...
            'Symbol', '', 'Labels', genotype_name_time_period);
    end
    hold on;
    scatter(repelem(1:length(genotype_select), cellfun(@length, frame_count_fly(b, :)))', horzcat(frame_count_fly{b, :}), 12); 
    ylabel(sprintf('# frames %s/fly', behav_name{b})); 
    xtickangle(45);
    hold off; 
    set(gcf, 'Renderer', 'painters'); 
    set(gcf, 'Position', get(0, 'Screensize'));
%     saveas(gcf, fullfile(output_folder, sprintf('%s-%s_frame_count.eps', dataset_name, behav_name{b})), 'epsc');
%     saveas(gcf, fullfile(output_folder, sprintf('%s-%s_frame_count.png', dataset_name, behav_name{b})));
end

if same_genotype_pairing
    bout_count_pair = cellfun(@(c) sum(reshape(c, 2,[])), bout_count_fly, 'UniformOutput', false);
    frame_count_pair = cellfun(@(c) sum(reshape(c, 2,[])), frame_count_fly, 'UniformOutput', false); 
    bout_count_pair_normalized = cell(size(bout_count_pair));
    frame_count_pair_normalized = cell(size(frame_count_pair)); 

    for g=1:length(genotype_select)
        for b=1:num_of_label
            n_bouts_per_fly = cellfun(@length, cluster_labels_cell_by_bout{g}); 
            bout_count_pair_normalized{b, g} = bout_count_pair{b, g}./sum(reshape(n_bouts_per_fly, 2, []));

            n_frames_per_fly = cellfun(@sum, cluster_bout_length_cell{g}); 
            frame_count_pair_normalized{b, g} = frame_count_pair{b, g}./sum(reshape(n_frames_per_fly, 2, []));
        end
    end
    
    for b=1:num_of_label
        figure()
        boxplot(horzcat(bout_count_pair{b, :}), repelem(1:length(genotype_select), cellfun(@length, bout_count_pair(b, :)))', ...
            'Symbol', '', 'Labels', genotype_name);
        hold on;
        scatter(repelem(1:length(genotype_select), cellfun(@length, bout_count_pair(b, :)))', horzcat(bout_count_pair{b, :}), 12); 
        ylabel(sprintf('# bouts %s/pair', behav_name{b})); 
        xtickangle(45);
        hold off; 
        set(gcf, 'Renderer', 'painters'); 
        set(gcf, 'Position', get(0, 'Screensize'));
%         saveas(gcf, fullfile(output_folder, sprintf('%s-%s_bout_count.eps', dataset_name, behav_name{b})), 'epsc');
%         saveas(gcf, fullfile(output_folder, sprintf('%s-%s_bout_count.png', dataset_name, behav_name{b})));

        figure()
        boxplot(horzcat(frame_count_pair{b, :}), repelem(1:length(genotype_select), cellfun(@length, frame_count_pair(b, :)))', ...
            'Symbol', '', 'Labels', genotype_name); 
        hold on;
        scatter(repelem(1:length(genotype_select), cellfun(@length, frame_count_pair(b, :)))', horzcat(frame_count_pair{b, :}), 12); 
        ylabel(sprintf('# frames %s/pair', behav_name{b})); 
        xtickangle(45);
        hold off; 
        set(gcf, 'Renderer', 'painters'); 
        set(gcf, 'Position', get(0, 'Screensize'));
%         saveas(gcf, fullfile(output_folder, sprintf('%s-%s_frame_count.eps', dataset_name, behav_name{b})), 'epsc');
%         saveas(gcf, fullfile(output_folder, sprintf('%s-%s_frame_count.png', dataset_name, behav_name{b})));
    end
    
    bout_count_fly = cell2struct(bout_count_fly, genotype_fieldname, 2);
    bout_count_fly_normalized = cell2struct(bout_count_fly_normalized, genotype_fieldname, 2);
    frame_count_fly = cell2struct(frame_count_fly, genotype_fieldname, 2);
    frame_count_fly_normalized = cell2struct(frame_count_fly_normalized, genotype_fieldname, 2);

    % save(fullfile(output_folder, sprintf('count_per_fly-%s.mat', dataset_name)), 'bout_count_fly', 'bout_count_fly_normalized', 'frame_count_fly', 'frame_count_fly_normalized');


    bout_count_pair = cell2struct(bout_count_pair, genotype_fieldname, 2);
    bout_count_pair_normalized = cell2struct(bout_count_pair_normalized, genotype_fieldname, 2);
    frame_count_pair = cell2struct(frame_count_pair, genotype_fieldname, 2);
    frame_count_pair_normalized = cell2struct(frame_count_pair_normalized, genotype_fieldname, 2);

    save(fullfile(output_folder, sprintf('count_per_pair-%s.mat', dataset_name)), 'bout_count_pair', 'bout_count_pair_normalized', 'frame_count_pair', 'frame_count_pair_normalized');
end


%% visualize transition matrix
text_label_x_offset = -0.2;
num_of_clusters_postproc = 5; 
active_clusters = 1:5; 

mc_count_cell = cell(size(genotype_select)); 
figure();
for g=1:length(genotype_select)
    if isempty(time_period)
        subplot(1, length(genotype_select), g); 
    else
        subplot(length(genotype_name), length(time_period), g); 
    end
    [markov_transition_by_fly, markov_transition_overall, markov_transition_reverse_overall] = ...
        analyze_markov_transition_helper(cluster_labels_cell_by_bout{g}, num_of_clusters_postproc, active_clusters);

    markov_transition_wo_normalize = sum(cat(3, markov_transition_by_fly{:}), 3);
    start_end_bout = cellfun(@(c) [c(1), c(end)], cluster_labels_cell_by_bout{g}, 'UniformOutput', false); 
    start_end_bout = vertcat(start_end_bout{:}); 
    for b=1:num_of_clusters_postproc
        assert(sum(markov_transition_wo_normalize(b, :)) - sum(markov_transition_wo_normalize(:, b)) == ...
            -diff(sum(start_end_bout == b, 1)))
    end

    imagesc(markov_transition_overall); 
    colorbar; 
    caxis([0, 1]);
    xticks(1:length(active_clusters));
    yticks(1:length(active_clusters));
    if isempty(time_period) 
        title(genotype_name{g}); 
    else
        [gt, gg] = ind2sub([length(time_period), length(genotype_name)], g); 
        title(sprintf('%s %s', genotype_name{gg}, time_period{gt})); 
    end
    ylabel('current behavior bout'); 
    xlabel('next behavior bout'); 

    text_xloc = repelem(1:length(active_clusters), length(active_clusters))+text_label_x_offset;
    text_yloc = repmat(1:length(active_clusters), 1, length(active_clusters));
    text_contents = strcat(cellstr(num2str(round(markov_transition_overall(:)*100, 1))), '%');
    text_contents = cellfun(@(x,y) {x;strcat('(', num2str(y), ')')}, text_contents, num2cell(markov_transition_wo_normalize(:)), 'UniformOutput', false);
    text(text_xloc, text_yloc, text_contents, 'FontSize', 8);

    set(gcf,'renderer','Painters');
    axis('square'); 
%     set(gcf, 'Position', [4, 45, 1236, 933]);
    set(gcf, 'Position', get(0, 'Screensize'));
    mc_count_cell{g} = markov_transition_wo_normalize; 
end

% saveas(gcf, fullfile(output_folder, sprintf('transition_matrix-%s-%d_%d.eps', dataset_name, genotype_select{1}, genotype_select{2})), 'epsc'); 
% saveas(gcf, fullfile(output_folder, sprintf('transition_matrix-%s-%d_%d.png', dataset_name, genotype_select{1}, genotype_select{2}))); 


% behav_category_struct = struct();
% for c=1:length(behav_name)-1
%     behav_category_struct.(matlab.lang.makeValidName(behav_name{c})) = c; 
% end
% figure();
% xloc = [0, 0.5, 0.5, 1.5, 1.5];
% yloc = [0, sqrt(3)/2, -sqrt(3)/2, sqrt(3)/2, -sqrt(3)/2];
% text_label_x_offset = 0.25.*[-1, -1, -1, 1/3, 1/3];
% for g=1:length(genotype_select)
%     subplot(1, 2, g);
%     [~, dp] = visualize_digraph(num2cell(active_clusters', 1), mc_count_cell{g}, 1, behav_category_struct, lines(length(behav_name)-1), 0.05, ...
%         xloc, yloc);
%     xticks([]);
%     yticks([]);
%     title(genotype_name{g}); 
%     dp.NodeLabel = behav_name(1:end-1);
% 
%     set(gcf,'renderer','Painters');
% end

%% visualize difference in transition matrix and test for statistical significance
if length(genotype_select) == 2
    text_label_x_offset = -0.2; 

    figure();
    markov_transition_diff = (mc_count_cell{2}./sum(mc_count_cell{2}, 2)) - (mc_count_cell{1}./sum(mc_count_cell{1}, 2));
    imagesc(markov_transition_diff);
    colormap(redbluecmap);
    caxis([-0.25, 0.25]);
    colorbar;
    xticks(1:length(active_clusters));
    yticks(1:length(active_clusters));
    title({genotype_name{2}, 'minus', genotype_name{1}}); 
    ylabel('current behavior bout'); 
    xlabel('next behavior bout'); 
 
    markov_transition_diff_synthetic = zeros([size(markov_transition_overall), n_repeat]);

%     n_fly_min = min(cellfun(@length, cluster_labels_cell_by_bout));
%     cluster_labels_cell_by_bout = cellfun(@(f) f(1:n_fly_min), cluster_labels_cell_by_bout, 'UniformOutput', false); 
%     cluster_bout_length_cell = cellfun(@(f) f(1:n_fly_min), cluster_bout_length_cell, 'UniformOutput', false);
%     cluster_labels_changepoint_cell = cellfun(@(f) f(1:n_fly_min), cluster_labels_changepoint_cell, 'UniformOutput', false); 

    for r=1:n_repeat
%         rand_idx = randperm(2*n_fly_min); 
        rand_idx = randperm(sum(cellfun(@length, cluster_labels_cell_by_bout))); 
        cluster_labels_cell_by_bout_all = horzcat(cluster_labels_cell_by_bout{:});
        cluster_labels_cell_by_bout_all = cluster_labels_cell_by_bout_all(rand_idx);
%         cluster_labels_cell_by_bout_permute = mat2cell(cluster_labels_cell_by_bout_all, 1, repmat(n_fly_min, 1, 2)); 
        cluster_labels_cell_by_bout_permute = mat2cell(cluster_labels_cell_by_bout_all, 1, cellfun(@length, cluster_labels_cell_by_bout)); 

        markov_transition_overall_permute_temp = zeros([size(markov_transition_overall), 2]);
        for g=1:length(genotype_select)
            [~, markov_transition_overall_permute_temp(:, :, g), ~] = ...
                analyze_markov_transition_helper(cluster_labels_cell_by_bout_permute{g}, 6, active_clusters);
        end
        markov_transition_diff_permute = markov_transition_overall_permute_temp(:, :, 2) - markov_transition_overall_permute_temp(:, :, 1);
        markov_transition_diff_synthetic(:, :, r) = markov_transition_diff_permute; 
    end
    markov_transition_diff_synthetic = sort(markov_transition_diff_synthetic, 3); 
    markov_transition_signif_thres_low = markov_transition_diff_synthetic(:, :, n_tail); 
    markov_transition_signif_thres_high = markov_transition_diff_synthetic(:, :, n_repeat-n_tail); 

    markov_transition_diff_signif = zeros(size(markov_transition_diff)); 
    markov_transition_diff_signif_mark = cell(size(markov_transition_diff_signif));
    for p=1:numel(markov_transition_diff)
        if markov_transition_diff(p) >= 0
            markov_transition_diff_signif(p) = markov_transition_diff(p) > markov_transition_signif_thres_high(p); 
        else
            markov_transition_diff_signif(p) = markov_transition_diff(p) < markov_transition_signif_thres_low(p); 
        end

        if markov_transition_diff_signif(p)
            markov_transition_diff_signif_mark{p} = '*'; 
        else
            markov_transition_diff_signif_mark{p} = ' '; 
        end

        if p == 9 || p == 14 || p == 17
            [r, c] = ind2sub(size(markov_transition_diff), p); 
            idx = find(markov_transition_diff(p) < squeeze(markov_transition_diff_synthetic(r, c, :)), 1); 
            if isempty(idx)
                idx = 10001; 
            end
            fprintf("%d-%d\n", p, 10001 - idx); 
        end
    end

    text_xloc = repelem(1:length(active_clusters), length(active_clusters))+text_label_x_offset;
    text_yloc = repmat(1:length(active_clusters), 1, length(active_clusters));
    text_contents = strcat(cellstr(num2str(round(markov_transition_diff(:)*100, 1))), '%');
    text_contents = cellfun(@(x,y) {x;y}, text_contents, markov_transition_diff_signif_mark(:), 'UniformOutput', false);
    text(text_xloc, text_yloc, text_contents, 'FontSize', 8);

    set(gcf,'renderer','Painters');
    % set(gcf, 'Position', get(0, 'Screensize'));
    axis('square'); 
    set(gcf, 'Position', [562, 503, 678, 475]); 
%     saveas(gcf, fullfile(output_folder, sprintf('transition_matrix_difference-%s-%d_%d.eps', dataset_name, genotype_select{1}, genotype_select{2})), 'epsc'); 
%     saveas(gcf, fullfile(output_folder, sprintf('transition_matrix_difference-%s-%d_%d.png', dataset_name, genotype_select{1}, genotype_select{2}))); 
end
