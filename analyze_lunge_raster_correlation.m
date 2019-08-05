close all;
% flymat_name = 'FLYMAT_CSMH_SM_HeisenbergChamber_Apr-May2019.mat';
% common_path = 'Z:\Kenichi\CSMH_SM_HeisenbergChamber_Apr-May2019';

flymat_name = 'FLYMAT_CG3385_KO_CS-bc-11_HeisenbergChamber_Apr2019.mat';
common_path = 'Z:\Kenichi\CG3385_KO_CS-bc-11_HeisenbergChamber_Apr-May2019';

load(fullfile(common_path, flymat_name));

max_lag_sec = 15;
fps = 60;
label_secs = linspace(-max_lag_sec, max_lag_sec, 2*10+1)';
smooth_window_width = fps; 
% Compute cross-correlation of lunge raster of two flies in a pair
% on the entire duration of experiment

crosscorr_roc_diff = zeros(length(flymatAll)/2, 1);
on_pad_time = zeros(length(flymatAll)/2, 2);
for i=1:2:length(flymatAll)
    lunge_raster_cell = cell(2, 1);
    
    % Compute auto-correlation for each fly's lunge raster
    for j=i:i+1
        rel_idx = 2 - mod(j, 2);
        % lunge_raster is a binary discrete sequence which is 1 at the
        % start of lunges
        % Alternatively, use flymatAll(j).L_binary hints on the
        % distribution of lunge durations
        lunge_raster_cell{rel_idx} = sparse(flymatAll(j).L_startsm, 1, 1, length(flymatAll(j).L_binary), 1);
        if isempty(lunge_raster_cell{rel_idx})
            lunge_raster_cell{rel_idx} = zeros(size(flymatAll(j).L_binary));
        end
        [autocorr_val, lags] = xcorr(full(lunge_raster_cell{rel_idx}), max_lag_sec*fps, 'biased');
        autocorr_val(ceil(length(autocorr_val)/2)) = 0;
        autocorr_val_sm = smoothdata(autocorr_val, 'gaussian', smooth_window_width); % Smooth correlation with Gaussian window 1sec wide
%         figure();
%         plot(lags, autocorr_val_sm);
%         title(sprintf('Autocorrelation of Fly %d Movie %s', flymatAll(j).fly, ...
%             strrep(flymatAll(j).movie{1}, '_', '-')));
%         xticks(label_secs.*fps);
%         xticklabels(cellstr(num2str(label_secs)));
%         xlabel('Lag (sec)');
    end
    
    % Compute cross-correlation between the two flies' lunge rasters
    [corr_val, lags] = xcorr(full(lunge_raster_cell{1}), full(lunge_raster_cell{2}), max_lag_sec*fps, 'biased');
    corr_val_sm = smoothdata(corr_val, 'gaussian', smooth_window_width);  % Smooth correlation with Gaussian window 1sec wide
%     figure();
%     plot(lags, corr_val_sm);
%     title(sprintf('Crosscorrelation between Fly %d and Fly %d of Movie %s', ...
%         flymatAll(i).fly, flymatAll(i+1).fly, strrep(flymatAll(j).movie{1}, '_', '-')));
%     xticks(label_secs.*fps);
%     xticklabels(cellstr(num2str(label_secs)));
%     xlabel('Lag (sec)');
    
    crosscorr_roc_diff(ceil(i/2)) = trapz(corr_val_sm(ceil(length(corr_val_sm)/2):end)) - trapz(corr_val_sm(1:floor(length(corr_val_sm/2))));
    
    % Identify food pad in the arena and compute on-pad time
    % Assume there is always four flies in a Heisenberg chamber
    if mod(i, 4) == 1
        movie_name_elems = strsplit(flymatAll(i).movie{1}, '_');
        movie_path = fullfile(common_path, movie_name_elems{1}, flymatAll(i).movie{1}, strcat(flymatAll(i).movie{1}, '.avi'));
        vin = VideoReader(movie_path);
        init_frame = readFrame(vin);
        init_frame = rgb2gray(init_frame);
        num_of_thres = 2; 
        thres = multithresh(init_frame, num_of_thres);
        init_frame_quant = imquantize(init_frame, thres);
        init_frame_binary = init_frame_quant == num_of_thres;
        s_elem = strel('square', 30);
        init_frame_binary = imclose(init_frame_binary, s_elem);
        conn_comps = bwconncomp(init_frame_binary);
        conn_comps_sz_list = cellfun(@length, conn_comps.PixelIdxList);
        [~, conn_comps_indices] = sort(conn_comps_sz_list, 'descend');
        food_pads_cell = conn_comps.PixelIdxList(conn_comps_indices(2:3));
        if max(food_pads_cell{1}) > max(food_pads_cell{2})
            food_pads_cell = food_pads_cell([2, 1]);
        end
        
        trx_path = fullfile(common_path, movie_name_elems{1}, flymatAll(i).movie{1}, flymatAll(i).movie{1}, strcat(flymatAll(i).movie{1}, '_JAABA'), 'trx.mat');
        load(trx_path);
    end
    
    food_pad_idx = ceil(abs(1.5 - mod(i, 4)));
    for j=i:i+1
        rel_idx = 2 - mod(j, 2);
        track_one_fly = [trx(flymatAll(j).fly).x, trx(flymatAll(j).fly).y];
        track_one_fly = floor(track_one_fly) + 1;
        track_one_fly_idx = sub2ind(size(init_frame), track_one_fly(:, 2), track_one_fly(:, 1));
        on_pad_mask = ismember(track_one_fly_idx, food_pads_cell{food_pad_idx});
        on_pad_time(ceil(i/2), rel_idx) = sum(on_pad_mask);
    end
end
%%
figure();
bar(on_pad_time, 'stacked');
title('Time spent on food pad for each pair of flies');
xlabel('Fly pair No.');
ylabel('Time on food pad (frames)');
legend({'Fly 1', 'Fly 2'});

figure();
bar(diff(on_pad_time, 1, 2));
title('Difference in time spent on food pad for each pair');
xlabel('Fly pair No.');
ylabel('Difference in time on food pad (frames)');

figure();
bar(crosscorr_roc_diff);
title('Difference in cross-correlation curve area-under-curve for each pair');
xlabel('Fly pair No.');
ylabel('Difference in area-under-curve');

