filename = ''; % Ephys dataname
ch_num = 16;
data = read_data(filename,ch_num);
% Filter Ephys data & extract bouts
[filtdata_r,filtdata_r_denoised] = filter_data(data.ch1); % right channel
[filtdata_l,filtdata_l_denoised] = filter_data(data.ch2); % left channel
bouts = extract_bouts(filtdata_l_denoised,filtdata_r_denoised);
%% Segment trials 
[trialData,idx] = segment_trials_MODI(data, bouts);
idx = round(idx,2);