load('Fig4g_data_HM_susI.mat');
load('Fig4g_data_HM_transI.mat');
%%
figure;
plotShadedErrorbar(HM_susI(mean(HM_susI(:,size(HM_susI,2)-5:size(HM_susI,2)-1),2)>0.3,1:size(HM_susI,2)-1),2.564,[255, 127, 14]/256,0.3);
plotShadedErrorbar(HM_transI(:,1:size(HM_susI,2)-1)+0.3,2.564,[44, 160, 44]/256,0.3);
ylabel('Normalized DFF');
xlabel('Time (s)');
set(gca,'FontSize',16);
xline(40/2.564);

legend({'','Sustained I','','Transient I',''});
legend boxoff
xlim([0 26]);





