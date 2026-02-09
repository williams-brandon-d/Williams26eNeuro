function IDs = removeBadLFP(IDs,info,dataSet)

        N = numel(IDs);

        remove_flag = false(N,1);

        for i = 1:N

            fileID = IDs{i};
            fileID_index = find_index(info,'ID',fileID);
    
            cell_info = info(fileID_index);
    
            % same cell num for dnqx before and after 
            cell_num = {};
            switch dataSet
                case 'Camk2'
                    switch cell_info.cell_type
                        case 'stellate'
                            cell_num = {'2','3','8'}; % 2 and 3 have low SNR, 8 has artefacts
                        case 'pyramidal'
                            cell_num = {'3','8','17'}; % 3 has strong 60 Hz, 8 has low SNR, 17 has artefacts
                    end

            end
            
            if any( ismember(cell_num,cell_info.cell_num) )
                remove_flag(i) = 1;
            end

        end

        IDs(remove_flag) = [];

end