load('ExtDataFig6a_integrator_distribution.mat');
%%
integ_cell_region_all = [];
for i = 1:length(integ_cell_region)
    integ_cell_region_all = [integ_cell_region_all,integ_cell_region{i}];
end

uniq_region_label = unique(integ_cell_region_all);
uniq_region_label = sort(uniq_region_label);
%
region_info = zeros(1,length(uniq_region_label));
region_ticks = {};
for i = 1:length(uniq_region_label)
    region_info(i) = length(find(integ_cell_region_all==uniq_region_label(i)));
    region_ticks{i} = region_label(uniq_region_label(i)+1).abbrev;
end
%
figure;
[~,I] = sort(region_info);
bar(region_info(flip(I)));
xticks(1:length(region_info));
xticklabels(region_ticks(flip(I)));
set(gca,'FontSize',16,'LineWidth',1.2);
ylabel('Cell number');
set(gcf,'Position',[130 1200 900 420]);
box off