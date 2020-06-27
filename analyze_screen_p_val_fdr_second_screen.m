excel_path = 'D:\xubo\NewTrainingFiles-VideoOnly\RNAi_screen_30fps';
excel_filename = '200625_2ndscreen_57genes-xubo_copy.xlsx';

% Read in 2nd screen data
rnai_data = readtable(fullfile(excel_path, excel_filename));
vdrc_line_name = rnai_data.Properties.VariableNames'; 
rnai_data = table2array(rnai_data)';

% Compile elav controls
elav_control_mask = cellfun(@(s) ~isempty(regexp(s, '.*(260b).*', 'tokens')), vdrc_line_name);
elav_control = rnai_data(elav_control_mask, :); 
elav_control = elav_control(~isnan(elav_control));

% Compile UAS controls
uas_control_mask = cellfun(@(s) ~isempty(regexp(s, '.*(cont).*', 'tokens')), vdrc_line_name);
uas_control = rnai_data(uas_control_mask, :); 
uas_control_cell = num2cell(uas_control, 2);
uas_control_cell = cellfun(@(d) d(~isnan(d)), uas_control_cell, 'UniformOutput', false);

control_mask = bitor(elav_control_mask, uas_control_mask);

% Compute rank-sum p-values and add to a new copy of the table
elav_p_val = nan(size(control_mask));
uas_p_val = nan(size(control_mask));
for i=1:length(control_mask)
    if control_mask(i)
        continue;
    else
        sample = rnai_data(i, :); 
        sample = sample(~isnan(sample));
        if isempty(sample)
            continue;
        else
            elav_p_val(i) = ranksum(sample, elav_control);
            uas_p_val(i) = ranksum(sample, uas_control_cell{(i-1)/2});
        end
    end
end

% Compile table of testing against elav control
rnai_mutant_against_elav_control = array2table(rnai_data, 'VariableNames', arrayfun(@(p) sprintf('Pair_%d', p), 1:size(rnai_data, 2), 'UniformOutput', false));
rnai_mutant_against_elav_control = addvars(rnai_mutant_against_elav_control, vdrc_line_name, 'Before', 'Pair_1');
rnai_mutant_against_elav_control = addvars(rnai_mutant_against_elav_control, elav_p_val, 'Before', 'Pair_1');
num_hypotheses = size(rnai_mutant_against_elav_control, 1) - sum(isnan(elav_p_val));

% Compute FDR and add to the table
rnai_mutant_against_elav_control = sortrows(rnai_mutant_against_elav_control, 'elav_p_val');
rnai_mutant_against_elav_control = rnai_mutant_against_elav_control(1:num_hypotheses, :);
elav_FDR = rnai_mutant_against_elav_control.elav_p_val .* size(rnai_mutant_against_elav_control, 1) ./ (1:size(rnai_mutant_against_elav_control, 1))';
rnai_mutant_against_elav_control = addvars(rnai_mutant_against_elav_control, elav_FDR, 'Before', 'Pair_1');

% Compute median and add to the table
lunge_count = table2array(rnai_mutant_against_elav_control(:, 4:end));
median_val = cellfun(@(counts) median(counts(~isnan(counts))), num2cell(lunge_count, 2));
rnai_mutant_against_elav_control = addvars(rnai_mutant_against_elav_control, median_val, 'Before', 'Pair_1');
num_pair_tested = cellfun(@(counts) sum(~isnan(counts)), num2cell(lunge_count, 2));
rnai_mutant_against_elav_control = addvars(rnai_mutant_against_elav_control, num_pair_tested, 'Before', 'Pair_1');


% Repeat for testing against UAS control
rnai_mutant_against_uas_control = array2table(rnai_data, 'VariableNames', arrayfun(@(p) sprintf('Pair_%d', p), 1:size(rnai_data, 2), 'UniformOutput', false));
rnai_mutant_against_uas_control = addvars(rnai_mutant_against_uas_control, vdrc_line_name, 'Before', 'Pair_1');
rnai_mutant_against_uas_control = addvars(rnai_mutant_against_uas_control, uas_p_val, 'Before', 'Pair_1');
num_hypotheses = size(rnai_mutant_against_uas_control, 1) - sum(isnan(uas_p_val));

rnai_mutant_against_uas_control = sortrows(rnai_mutant_against_uas_control, 'uas_p_val');
rnai_mutant_against_uas_control = rnai_mutant_against_uas_control(1:num_hypotheses, :);
uas_FDR = rnai_mutant_against_uas_control.uas_p_val .* size(rnai_mutant_against_uas_control, 1) ./ (1:size(rnai_mutant_against_uas_control, 1))';
rnai_mutant_against_uas_control = addvars(rnai_mutant_against_uas_control, uas_FDR, 'Before', 'Pair_1');

lunge_count = table2array(rnai_mutant_against_uas_control(:, 4:end));
median_val = cellfun(@(counts) median(counts(~isnan(counts))), num2cell(lunge_count, 2));
rnai_mutant_against_uas_control = addvars(rnai_mutant_against_uas_control, median_val, 'Before', 'Pair_1');
num_pair_tested = cellfun(@(counts) sum(~isnan(counts)), num2cell(lunge_count, 2));
rnai_mutant_against_uas_control = addvars(rnai_mutant_against_uas_control, num_pair_tested, 'Before', 'Pair_1');

writetable(rnai_mutant_against_elav_control, excel_filename, 'Sheet', 'FDR_mutant_against_elav_control');
writetable(rnai_mutant_against_uas_control, excel_filename, 'Sheet', 'FDR_mutant_against_uas_control');

% Combine two sets of test and repeat 
rnai_mutant_against_elav_control.Properties.VariableNames{'elav_p_val'} = 'p_val';
rnai_mutant_against_elav_control = removevars(rnai_mutant_against_elav_control, 'elav_FDR');
control_line = repmat({vdrc_line_name{1}}, height(rnai_mutant_against_elav_control), 1);
rnai_mutant_against_elav_control = addvars(rnai_mutant_against_elav_control, control_line, 'After', 'vdrc_line_name');
rnai_mutant_against_uas_control.Properties.VariableNames{'uas_p_val'} = 'p_val';
rnai_mutant_against_uas_control = removevars(rnai_mutant_against_uas_control, 'uas_FDR');
control_line = cellfun(@(s) strrep(s, 'elav', 'cont'), rnai_mutant_against_uas_control.vdrc_line_name, 'UniformOutput', false);
rnai_mutant_against_uas_control = addvars(rnai_mutant_against_uas_control, control_line, 'After', 'vdrc_line_name');

rnai_mutant_against_control_combined = [rnai_mutant_against_elav_control; rnai_mutant_against_uas_control];
num_hypotheses = size(rnai_mutant_against_control_combined, 1);

rnai_mutant_against_control_combined = sortrows(rnai_mutant_against_control_combined, 'p_val');
combined_FDR = rnai_mutant_against_control_combined.p_val .* size(rnai_mutant_against_control_combined, 1) ./ (1:size(rnai_mutant_against_control_combined, 1))';
rnai_mutant_against_control_combined = addvars(rnai_mutant_against_control_combined, combined_FDR, 'After', 'p_val');

% Save sheets
writetable(rnai_mutant_against_control_combined, excel_filename, 'Sheet', 'FDR_against_control_combined');
