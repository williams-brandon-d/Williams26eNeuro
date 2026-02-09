function filteredData = filterData(data,Fs,bandpass,order)
    hpf_fc = bandpass(1);
    lpf_fc = bandpass(2);

    if hpf_fc > 0
        [b_hpf,a_hpf] = butter(order,hpf_fc/(Fs/2),'high'); % hpf coefficients
        hpf_data = filtfilt(b_hpf,a_hpf,data); % hpf data
    else
        hpf_data = data;
    end
    
    if lpf_fc > 0
        [b_lpf,a_lpf] = butter(order,lpf_fc/(Fs/2),'low');
        lpf_data = filtfilt(b_lpf,a_lpf,hpf_data);
    else
        lpf_data = hpf_data;
    end
    
    filteredData = lpf_data;
end