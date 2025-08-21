load('Fig2c_data_nanparam_kstest2_p_value.mat');
%%
frame_rate = 2.564;
rank_size = mean(-1*log10(nonparam_p_value),2);
[~,diff_ind] = sort(rank_size);
diff_ind = setdiff(diff_ind,find(isnan(rank_size)),'stable');
diff_ind = flip(diff_ind);
figure;
HM = -1*log10(nonparam_p_value(diff_ind,:));
% HM(HM>1.3) = 2;
h = imagesc(HM);
CM = othercolor('RdPu9');
colormap(CM);
clim([-1*log10(0.05) 8]);


for i = 1:length(diff_ind)
    YLabels{i} = region_label(diff_ind(i)+1).abbrev;
end
yticks(1:length(diff_ind));
yticklabels(YLabels);
xticks([5 10]*frame_rate+0.5);
xticklabels({'5','10'})
line([[5 10]*frame_rate+0.5;[5 10]*frame_rate+0.5],[length(diff_ind),length(diff_ind)+2],'Color','k');
set(gca,'TickLength',[0 0]);
row = 2:37;
yline([row-.5, row+.5], 'w-'); % see footnotes [1,2]

set(gca,'FontSize',12);
xlabel('Time after obstacle disappears (s)');
colorbar;