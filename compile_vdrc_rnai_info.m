root_folder = ('D:\xubo\NewTrainingFiles-VideoOnly\RNAi_screen_30fps');
screen_excel = '200317_SupplementaryInformation-1.xlsx';
screen_excel_new = '200317_SupplementaryInformation-1-id_complete-option_4.xlsx';
sheet_names = {'FDR_old', 'FDR_new'};
complete_vdrc_info_excel = 'REPORT_VdrcCatalogueKK.xls';

for i=1:length(sheet_names)
    screen_table = readtable(fullfile(root_folder, screen_excel), 'Sheet', sheet_names{i});
    vdrc_info_table = readtable(fullfile(root_folder, complete_vdrc_info_excel));
    cg_number = cell(size(screen_table, 1), 1);
    construct_id = zeros(size(screen_table, 1), 1);
    synonyms = cell(size(screen_table, 1), 1);
    for j=1:size(screen_table, 1)
        vdrc_id = screen_table.VDRCID(j);
        vdrc_mask = vdrc_info_table.vdrcId == vdrc_id; 
        
        assert(sum(vdrc_mask) ~= 0, 'No corresponding VDRC line found');
        assert(sum(vdrc_mask) == 1, 'More than one corresponding VRDC line found');
        
        cg_number{j} = vdrc_info_table.cgNumber{vdrc_mask};
        construct_id(j) = vdrc_info_table.constructId(vdrc_mask);
        synonyms{j} = vdrc_info_table.synonyms{vdrc_mask};
    end
    
    if strcmp(sheet_names{i}, 'FDR_old')
        screen_table = addvars(screen_table, cg_number, 'Before', 'VDRCID', 'NewVariableNames', {'CGNumber'});
    elseif strcmp(sheet_names{i}, 'FDR_new')
        try 
            cg_number_match_mask = cellfun(@(a, b) strcmp(a, b), cg_number, screen_table.CGNumber);
            assert(all(cg_number_match_mask));
        catch
            vdrc_synonyms = vdrc_info_table.synonyms;
            cg_number_mismatch_idx = find(~cg_number_match_mask);
            for k=1:length(cg_number_mismatch_idx)
                vdrc_id = screen_table.VDRCID(cg_number_mismatch_idx(k));
                vdrc_mask = vdrc_info_table.vdrcId == vdrc_id; 
                assert(contains(vdrc_synonyms{vdrc_mask}, screen_table.CGNumber{cg_number_mismatch_idx(k)}), ...
                    'CG number %s and %s does not match and %s not found in synonyms');
            end
            % Option 1: if all mismatching CG numbers found in synonyms, change all CG numbers to those used by VDRC official info sheet
            screen_table.CGNumber = cg_number; 
        end
    end
    screen_table = addvars(screen_table, construct_id, 'After', 'VDRCID', 'NewVariableNames', {'construct_ID'});
    % Option 3: also include synonym column 
%     screen_table = addvars(screen_table, synonyms, 'After', 'construct_ID', 'NewVariableNames', {'synonym'});

    writetable(screen_table, fullfile(root_folder, screen_excel_new), 'Sheet', sheet_names{i});
    
    % Option 4: include synonym as separate sheet
    synonyms_table = array2table([cg_number, synonyms], 'VariableNames', {'CGNumber', 'synonyms'});
    writetable(synonyms_table, fullfile(root_folder, screen_excel_new), 'Sheet', sprintf('synonyms_%s', sheet_names{i}));
end
