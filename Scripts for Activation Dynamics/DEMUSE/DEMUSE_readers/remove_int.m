function signal_f = remove_int(signal,fsamp,cut_high,cut_low)

[b,a] = butter(2,[cut_high/(fsamp/2) cut_low/(fsamp/2)]);
signal_f = filtfilt(b,a,signal);
