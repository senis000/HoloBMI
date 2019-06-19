function fa_estParams2sot_neuron_direct_indirect...
    (estParams, direct_bool, direct_idxs, indirect_bool, indirect_idxs)
%4.12.17


%%



%%
            if(direct_bool && indirect_bool)
                num_direct_i = ...
                    length(chr2_data_mat_over_animal_session(animal_idx, session_idx).direct_idxs);
                num_indirect_i = ...
                    length(chr2_data_mat_over_animal_session(animal_idx, session_idx).indirect_idxs);

                num_direct_analyze_i = ...
                    sum(analyze_idxs_i(1:num_direct_i));
                num_indirect_analyze_i = ...
                    sum(analyze_idxs_i(num_direct_i+(1:num_indirect_i)));

                direct_sel      = 1:num_direct_analyze_i;
                indirect_sel    = num_direct_analyze_i+(1:num_indirect_analyze_i);
                
%             %For Debug:
%             disp(['num_total_analyze_i: ' num2str(num_total_analyze_i)]);
%             disp(['num_direct_analyze_i: ' num2str(num_direct_analyze_i)]);
%             disp(['num_indirect_analyze_i: ' num2str(num_indirect_analyze_i)]);
%             disp(['L size: ' num2str(size(estParams_i.L))]);
%             input('enter');                 
            
                %direct only:
                sel_idxs    = ...
                    direct_sel;
                [direct_result] = fa2sot_subset_of_channels(estParams_i, sel_idxs);
                %ASSIGN
                fa_results(animal_idx, session_idx).direct = ...
                    direct_result;
                SOT.direct(animal_idx, session_idx) = ...
                    direct_result.SOT_over_neuron_mean;

                %indirect only:
                sel_idxs    = ...
                    indirect_sel;
                [indirect_result] = fa2sot_subset_of_channels(estParams_i, sel_idxs);
                %ASSIGN
                fa_results(animal_idx, session_idx).indirect = ...
                    indirect_result;        
                SOT.indirect(animal_idx, session_idx) = ...
                    indirect_result.SOT_over_neuron_mean;
            
                %all:
                sel_idxs    = ...
                    [direct_sel indirect_sel]; 
                [all_result] = fa2sot_subset_of_channels(estParams_i, sel_idxs);
                %ASSIGN
                fa_results(animal_idx, session_idx).all = ...
                    all_result;        
                SOT.all(animal_idx, session_idx) = ...
                    all_result.SOT_over_neuron_mean;
            end
            
            if(indirect_bool)
                sel_idxs = 1:num_total_analyze_i;
                [indirect_result] = fa2sot_subset_of_channels(estParams_i, sel_idxs);
                %ASSIGN
                fa_results(animal_idx, session_idx).indirect = ...
                    indirect_result;        
                SOT.indirect(animal_idx, session_idx) = ...
                    indirect_result.SOT_over_neuron_mean;
            end