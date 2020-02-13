close all;
behavior_shorthands = {'L'};
max_lag_sec = 60;
fps = 60;
label_secs = linspace(-max_lag_sec, max_lag_sec, 2*10+1)';
smooth_window_width = fps/10; 
exp_length = 30;

exp_struct = struct('flymat_name', '', 'common_path', '', 'genotypes', {}); 

exp_struct(1).flymat_name = 'FLYMAT_CSMH_SM_HeisenbergChamber_Apr-May2019.mat';
exp_struct(1).common_path = 'Z:\Kenichi\CSMH_SM_HeisenbergChamber_Apr-May2019';
exp_struct(1).genotypes = {'CSMH\_SM', 'CSMH\_SM'};

exp_struct(2).flymat_name = 'FLYMAT_CSMH_SM_vsWT_GM_HeisenbergChamber_May2019.mat';
exp_struct(2).common_path = 'Z:\Kenichi\CSMH_SM_vsWT_GM_HeisenbergChamber_May2019';
exp_struct(2).genotypes = {'CSMH\_SM', 'WT\_GM'};

exp_struct(3).flymat_name = 'FLYMAT_WingThreat_bouts_ThirtyMinAfterTwoFliesEntered.mat';
exp_struct(3).common_path = 'Z:\Kenichi\CG3385_KO_CS-bc-11_GM_vsWT_GM_HeisenbergChamber_May2019';
exp_struct(3).genotypes = {'nervy\_mut', 'WT\_GM'};

autocorr_val_all_exp = cell(size(exp_struct));
dominant_fly_mask_all_exp = cell(size(exp_struct));

for m=1:length(exp_struct) 
    load(fullfile(exp_struct(m).common_path, exp_struct(m).flymat_name));
    
    autocorr_val_all = cell(length(flymatAll), length(behavior_shorthands));
    for i=1:2:length(flymatAll)
        behav_raster_cell = cell(2, length(behavior_shorthands));

        % Compute auto-correlation for each fly's lunge raster
        for j=i:i+1
            rel_idx = 2 - mod(j, 2);
            % lunge_raster is a binary discrete sequence which is 1 at the
            % start of lunges
            % Alternatively, use flymatAll(j).L_binary hints on the
            % distribution of lunge durations
            for k=1:length(behavior_shorthands)
                startsm_field = strjoin({behavior_shorthands{k}, 'startsm'}, '_');
                binary_field = strjoin({behavior_shorthands{k}, 'binary'}, '_');
                behav_raster_cell{rel_idx, k} = sparse(flymatAll(j).(startsm_field), 1, 1, length(flymatAll(j).(binary_field)), 1);
                if isempty(behav_raster_cell{rel_idx})
                    behav_raster_cell{rel_idx, k} = zeros(size(flymatAll(j).(binary_field)));
                end
%                 behav_raster_cell = cellfun(@(raster) normalize(raster, 'center'), behav_raster_cell, 'UniformOutput', false);

                [autocorr_val, lags] = xcorr(full(behav_raster_cell{rel_idx, k}), max_lag_sec*fps, 'coeff');
                if all(isnan(autocorr_val))
                    autocorr_val = zeros(size(autocorr_val));
                end
                autocorr_val(ceil(length(autocorr_val)/2)) = 0;
                autocorr_val_sm = smoothdata(autocorr_val, 'gaussian', smooth_window_width); % Smooth correlation with Gaussian window 1sec wide
%                 autocorr_val_sm = autocorr_val; 
                autocorr_val_all{i+rel_idx-1, k} = autocorr_val_sm; 
            end  
        end
    end
    
    lunge_count = arrayfun(@(s) length(s.L_start), flymatAll);
    lunge_count = reshape(lunge_count, 2, []);
    [~, dominant_fly_rel_idx] = max(lunge_count);
    dominant_fly_abs_idx = sub2ind(size(lunge_count), dominant_fly_rel_idx, 1:size(lunge_count, 2));
    dominant_fly_mask = false(size(flymatAll));
    dominant_fly_mask(dominant_fly_abs_idx) = true; 
    
    autocorr_val_all_exp{m} = autocorr_val_all;
    dominant_fly_mask_all_exp{m} = dominant_fly_mask; 
end
%%
figure();
hold on; 
legends = {}; 
ticks = 1:(max_lag_sec*fps+1); 
for m=1:length(exp_struct)
    autocorr_val_all = autocorr_val_all_exp{m};
    dominant_fly_mask = dominant_fly_mask_all_exp{m};
    
    if length(unique(exp_struct(m).genotypes)) > 1
        for i=1:2
            mean_autocorr = mean(horzcat(autocorr_val_all{i:2:length(autocorr_val_all)}), 2); 
            plot(mean_autocorr(ceil(length(mean_autocorr)/2):end));
%             f_obj = fit(ticks', mean_autocorr, 'poly1');
%             line(ticks-1, f_obj.p1.*ticks+f_obj.p2, 'HandleVisibility','off');
            legends = [legends, sprintf('%s in %s-v-%s', ...
                exp_struct(m).genotypes{i}, exp_struct(m).genotypes{1}, exp_struct(m).genotypes{2})];
        end
    else
        mean_autocorr = mean(horzcat(autocorr_val_all{dominant_fly_mask}), 2); 
        plot(mean_autocorr(ceil(length(mean_autocorr)/2):end));
%         f_obj = fit(ticks', mean_autocorr, 'poly1');
%         line(ticks-1, f_obj.p1.*ticks+f_obj.p2, 'HandleVisibility','off');
        legends = [legends, sprintf('dominant %s in %s-v-%s', ...
            exp_struct(m).genotypes{1}, exp_struct(m).genotypes{1}, exp_struct(m).genotypes{2})];
        
        mean_autocorr = mean(horzcat(autocorr_val_all{~dominant_fly_mask}), 2); 
        plot(mean_autocorr(ceil(length(mean_autocorr)/2):end));
%         f_obj = fit(ticks', mean_autocorr, 'poly1');
%         line(ticks-1, f_obj.p1.*ticks+f_obj.p2, 'HandleVisibility','off');
        legends = [legends, sprintf('subordinate %s in %s-v-%s', ...
            exp_struct(m).genotypes{1}, exp_struct(m).genotypes{1}, exp_struct(m).genotypes{2})];
    end
end

xlim([0, max_lag_sec*fps]);
xticks(0:2*fps:max_lag_sec*fps);
xticklabels(0:2:max_lag_sec);
xlabel('Time lag between two lunges (in s)');
ylabel('Mean smoothed autocorrelation value (a.u.)');

ax = gca; 
ax.ColorOrderIndex = 1; 
for m=1:length(exp_struct)
    ax = gca;
    plot(ticks, ones(size(ticks))./sqrt((length(autocorr_val_all_exp{m})/2)*exp_length*60*fps));
    ax.ColorOrderIndex = ax.ColorOrderIndex + 1;
    legends = [legends, sprintf('Standard error in %s-v-%s', exp_struct(m).genotypes{1}, exp_struct(m).genotypes{2})];
end

legend(legends);

