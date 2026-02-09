function file = detectArtifacts(file)
% alternatively could add theta and gamma components to make threshold

switch file.stimType
    case 'theta'
        % detect artifacts for all cycles in data 
        maxCycles = numel(file.cycle_start_index);
        
        for i = 1:file.nDataChannels
          fieldName = file.dataChannels{i};
%           file.(fieldName).spikes = false(file.nCycles,file.nTrials);
          file.(fieldName).spikes = false(maxCycles,file.nTrials);
        
            % skip cycles with spikes
            for iTrial = 1:file.nTrials 
                trial = file.trials(iTrial);
                trial_data = file.(fieldName).raw_data(:,trial);
        
                for iCycle = 1:maxCycles
%                     cycle = file.cycles(iCycle);
                    cycle = iCycle;
                    window_indices = (0:file.cycle_length) + file.cycle_start_index(cycle);
                    window_data = trial_data(window_indices);
        
                    if (range(window_data) > file.(fieldName).range_threshold)
                      file.(fieldName).spikes(cycle,trial) = 1;
                      fprintf('Spike Detected: Cycle %d, Trial %d\n',cycle,trial)
                    end
        
                end
        
            end
        
        end

    case 'pulse'

        for i = 1:file.nDataChannels
          fieldName = file.dataChannels{i};
          file.(fieldName).spikes = false(file.nTrials,1);
        
            for iTrial = 1:file.nTrials 
                trial = file.trials(iTrial);
                trial_data = file.(fieldName).raw_data(:,trial);
        
                window_indices = (0:file.pulse_length) + file.pulse_start_index;
                window_data = trial_data(window_indices);
    
                if (range(window_data) > file.(fieldName).range_threshold)
                  file.(fieldName).spikes(iTrial) = 1;
                  fprintf('Spike Detected: Trial %d\n',trial)
                end
        
            end
        
        end

end


end