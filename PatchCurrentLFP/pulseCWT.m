function [x,y,z] = pulseCWT(data,fb,time,type)

nTrials = size(data,2);

    switch type
        case 'average data'
            mean_filt_data = mean(data,2); % average all traces
            [wt,freq] = cwt(mean_filt_data,'FilterBank',fb);
            z = abs( wt ).^2;

        case 'average scalogram' % not fixed yet
            for iTrial = 1:nTrials
                trial_filt_data = data(:,iTrial);
            
                [wt,freq] = cwt(trial_filt_data,'FilterBank',fb);
                z_trial = abs( wt ).^2;
                
                if iTrial == 1
                    z_sum = z_trial;
                else
                    z_sum = z_sum + z_trial;
                end
            end
            
            z = z_sum/nTrials;
    end

    [x,y] = meshgrid(time,freq);


end