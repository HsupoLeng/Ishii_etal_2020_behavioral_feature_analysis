excel_path = 'D:\xubo\NewTrainingFiles-VideoOnly\RNAi_screen_30fps';
excel_filename = '200314_1stScreen-30fpsMovies_LungeNumbers_XL-KI_identical_lines_combined.xlsx';
old_excel_path = 'D:\xubo\NewTrainingFiles-VideoOnly\RNAi_screen_30fps\Kenta_RNAi_screen';
old_excel_filename = 'RNAiScreen1st-Old-All_identical_lines_combined.xlsx';
old_new_combine = 2; % For test on only old data, set 0; for test on only new data, set 1; for combined test, set 2

% Read in the new screening data
rnai_data = readtable(fullfile(excel_path, excel_filename));
lunge_count_cols_ori = 5:28;
rnai_data = rnai_data(:, 1:lunge_count_cols_ori(end));

% Compile control in the new data
control_mask = cellfun(@(s) ~isempty(regexp(s, '.*(260B).*', 'tokens')), rnai_data.Genotype);
control = rnai_data(control_mask, lunge_count_cols_ori); % Assume lunge counts for wells are in columns 5 to 16 
control = table2array(control);
control = control(~isnan(control));
control_vec = control(:);

% Read in old screening data and combine tables
rnai_data_old = readtable(fullfile(old_excel_path, old_excel_filename));
control_old = table2array(rnai_data_old(1:967, 2));
rnai_data_old = rnai_data_old(:, [1, 3:end]);
rnai_data_old = rnai_data_old(:, [true, ~all(isnan(table2array(rnai_data_old(:, 2:end))))]);
old_genotype = rnai_data_old.Properties.VariableNames(2:end)';
old_genotype = cellfun(@(s) strsplit(s, '_'), old_genotype, 'UniformOutput', false);
old_genotype = cellfun(@(c) erase(c{1}, ' '), old_genotype, 'UniformOutput', false);
old_genotype = cellfun(@(c) erase(c, 'x'), old_genotype, 'UniformOutput', false);
rnai_data_old_transp_genotype = array2table(old_genotype, 'VariableNames', {'Genotype'});

old_transp_count = table2array(rnai_data_old(:, 2:end))';
rnai_data_old_transp_count = array2table(old_transp_count, 'VariableNames', [rnai_data_old.Genotype(:)]);
rnai_data_old_transp_count = rnai_data_old_transp_count(:, ~all(isnan(old_transp_count)));
rnai_data_old_transp = [rnai_data_old_transp_genotype, rnai_data_old_transp_count];

new_varnames_in_rnai_data = rnai_data_old_transp.Properties.VariableNames(~ismember(rnai_data_old_transp.Properties.VariableNames, rnai_data.Properties.VariableNames));
for i=1:length(new_varnames_in_rnai_data)
    rnai_data.(new_varnames_in_rnai_data{i}) = nan(size(rnai_data, 1), 1);
end

new_varnames_in_rnai_data_old_transp = rnai_data.Properties.VariableNames(~ismember(rnai_data.Properties.VariableNames, rnai_data_old_transp.Properties.VariableNames));
for i=1:length(new_varnames_in_rnai_data_old_transp)
    rnai_data_old_transp = addvars(rnai_data_old_transp, repmat({''}, size(rnai_data_old_transp, 1), 1), 'After', i, 'NewVariableNames', new_varnames_in_rnai_data_old_transp{i});
end

if old_new_combine == 2
    rnai_data_all = [rnai_data; rnai_data_old_transp];
    control_all = [control_vec; control_old];
    control_mask = [control_mask; false(size(rnai_data_all, 1) - length(control_mask), 1)];
    lunge_count_cols = 5:size(rnai_data_all, 2);
elseif old_new_combine == 1
    rnai_data_all = rnai_data;
    control_all = control_vec;
    control_mask = control_mask; 
    lunge_count_cols = lunge_count_cols_ori; 
elseif old_new_combine == 0
    rnai_data_all = rnai_data_old_transp; 
    control_all = control_old;
    control_mask = false(size(rnai_data_all, 1), 1);
    lunge_count_cols = 5:size(rnai_data_all, 2);
end

% Compute rank-sum p-values and add to a new copy of the table
p_val = nan(size(control_mask));
for i=1:length(control_mask)
    if control_mask(i)
        continue;
    else
        sample = table2array(rnai_data_all(i, lunge_count_cols)); 
        sample = sample(~isnan(sample));
        if isempty(sample)
            continue;
        else
            p_val(i) = ranksum(sample, control_all);
        end
    end
end

rnai_data_all_new = addvars(rnai_data_all, p_val, 'After', 'Note');
num_hypotheses = size(rnai_data_all_new, 1) - sum(isnan(p_val));

% Compute FDR and add to new table
rnai_data_all_new = sortrows(rnai_data_all_new, 'p_val');
rnai_data_all_new = rnai_data_all_new(1:num_hypotheses, :);
FDR = rnai_data_all_new.p_val .* size(rnai_data_all_new, 1) ./ (1:size(rnai_data_all_new, 1))';
rnai_data_all_new = addvars(rnai_data_all_new, FDR, 'After', 'p_val');

lunge_count_cols = lunge_count_cols + 2;
% Compute median and and number of pairs tested, and add to new table
lunge_count = table2array(rnai_data_all_new(:, lunge_count_cols));
median_val = cellfun(@(counts) median(counts(~isnan(counts))), num2cell(lunge_count, 2));
rnai_data_all_new = addvars(rnai_data_all_new, median_val, 'After', 'FDR');
num_pair_tested = cellfun(@(counts) sum(~isnan(counts)), num2cell(lunge_count, 2));
rnai_data_all_new = addvars(rnai_data_all_new, num_pair_tested, 'After', 'median_val');


if old_new_combine == 2
    writetable(rnai_data_all_new, excel_filename, 'Sheet', 'FDR_all');
elseif old_new_combine == 1
    writetable(rnai_data_all_new, excel_filename, 'Sheet', 'FDR_new');
elseif old_new_combine == 0
    writetable(rnai_data_all_new, excel_filename, 'Sheet', 'FDR_old');
end

