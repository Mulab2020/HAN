load('ExtDataFig2c_ephys_recording.mat');
%%
figure;
x = (1:length(filtdata_l_denoised(1:10:end)))/600;
plot(x,filtdata_l_denoised(1:10:end));
hold on
plot(x,-filtdata_r_denoised(1:10:end));
hold off
set(gcf,"Position",[0 1109 2298 420]);