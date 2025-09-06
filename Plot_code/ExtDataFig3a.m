load('ExtDataFig3a_brainregion_PC1_resp.mat');
%%
frame_rate = 2.7;
colorgradient = othercolor('PiYG4');
colors = [colorgradient(1,:);colorgradient(256,:)];
figure;
% [ha, pos] = tight_subplot(7,8,[.02 .01],[.02 .02],[.02 .02]);
for i = 1:53
    
%     axes(har(i));
    subplot(7,8,i);
    title(region_label(i+1).abbrev);
    if isempty(PC_1_intvl_resp{i})
        continue;
    else
        hold on;
        PC_1_l_intvl = PC_1_intvl_resp{i}{1};
        PC_1_r_intvl = PC_1_intvl_resp{i}{2};
        plotShadedErrorbar(PC_1_l_intvl(:,1:end-1),frame_rate,colors(1,:),0.2);
        plotShadedErrorbar(PC_1_r_intvl(:,1:end-1),frame_rate,colors(2,:),0.2);
        set(gca,'Color','none', 'FontName', 'Calibri');
        xlim([0 (size(PC_1_l_intvl,2)-3)/frame_rate]);
        xticks([0 5 10]);
        box off
        hold off;

    end
end
set(findall(gcf, 'Type', 'text'), 'FontName', 'Calibri');
set(gcf,'position',[300   300   1200   1000]);