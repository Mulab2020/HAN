load('Fig3g_data_early_late_distribution');
%%
for region = 17
    if ~isnan(type1_resp_regions{region}(1,1))
        close all;
        type1_resp = type1_resp_regions{region};
        type2_resp = type2_resp_regions{region};
        type3_resp = type3_resp_regions{region};
        type5_resp = type5_resp_regions{region};
        type6_resp = type6_resp_regions{region};
        type7_resp = type7_resp_regions{region};
        % figure;
        % hold on;
        % errorbar(1,mean(type3_resp(:,sframe)),sem(type3_resp(:,sframe)));
        % errorbar(2,mean(type2_resp(:,sframe)),sem(type2_resp(:,sframe)));
        % errorbar(3,mean(type1_resp(:,sframe)),sem(type1_resp(:,sframe)));
        % errorbar(4,mean(type7_resp(:,sframe)),sem(type7_resp(:,sframe)));
        % errorbar(5,mean(type6_resp(:,sframe)),sem(type5_resp(:,sframe)));
        % errorbar(6,mean(type5_resp(:,sframe)),sem(type6_resp(:,sframe)));

        % -------------------------------------
        colorgradient = othercolor('PiYG4');

        colors = [colorgradient(1,:);colorgradient(52,:);colorgradient(103,:);colorgradient(154,:);colorgradient(205,:);colorgradient(256,:)];

        freq = (14077/(80*60));
        sframe = ceil(2*freq); % trial 17s (late delay)
        % region = 9;
        type1_sframe_resp = -type1_resp(:,sframe); type1_sframe_resp = type1_sframe_resp(:);
        type2_sframe_resp = -type2_resp(:,sframe); type2_sframe_resp = type2_sframe_resp(:);
        type3_sframe_resp = -type3_resp(:,sframe); type3_sframe_resp = type3_sframe_resp(:);
        type5_sframe_resp = -type5_resp(:,sframe); type5_sframe_resp = type5_sframe_resp(:);
        type6_sframe_resp = -type6_resp(:,sframe); type6_sframe_resp = type6_sframe_resp(:);
        type7_sframe_resp = -type7_resp(:,sframe); type7_sframe_resp = type7_sframe_resp(:);
%         L_mean_resp = mean([type1_sframe_resp;type2_sframe_resp;type3_sframe_resp],'omitnan');
%         R_mean_resp = mean([type5_sframe_resp;type6_sframe_resp;type7_sframe_resp],'omitnan');
        L_R_resp = [type1_sframe_resp;type2_sframe_resp;type3_sframe_resp;type5_sframe_resp;type6_sframe_resp;type7_sframe_resp];
        bias = (prctile(L_R_resp,3)+prctile(L_R_resp,98))/2;
        scaling = abs((prctile(L_R_resp,98)-prctile(L_R_resp,3))/2);
        all_resp = ([type1_sframe_resp;type2_sframe_resp;type3_sframe_resp;type1_sframe_resp;type6_sframe_resp;type7_sframe_resp]-bias)/scaling;
        all_resp_cell = {(type3_sframe_resp-bias)/scaling;(type2_sframe_resp-bias)/scaling;(type1_sframe_resp-bias)/scaling;(type7_sframe_resp-bias)/scaling;(type6_sframe_resp-bias)/scaling;(type5_sframe_resp-bias)/scaling};
        figure;
        % ax = axes;
        subplot(1,3,1);
        hold on;
        
        range = linspace(min(all_resp),max(all_resp),15);
        range = [2*range(1)-range(2),range,2*range(end)-range(end-1)];
        % range = min(all_resp):0.0015:max(all_resp);
        
        for i = 1:6
            EBdata = all_resp_cell{i};
            B = 5000;
            histfun = @(x) histcounts(x, range)/length(x);
            bootstrap_proportion = bootstrp(B, histfun, EBdata);
            %     sem_counts = std(bootstrap_counts);
            %     sem_proportions = sem_counts / length(EBdata);

            plotShadedErrorbar_flip_givenY((range(1:end-1)+range(2:end))/2,bootstrap_proportion,colors(i,:),0.3);
            yline(mean(EBdata,'omitnan'),'Color',colors(i,:),'Linestyle','--','Linewidth',1.5);

        end
        yticks([-1 0 1]);
        ylim([-1 1]);
        
        set(gca,'FontSize',16,'LineWidth',0.5);
        set(findall(gcf, 'Type', 'text'), 'FontName', 'Calibri'); % 所有文本

%         
        

        % ------------------------------------------
        sframe = ceil(17*freq); % trial 17s (late delay)
        type1_sframe_resp = -type1_resp(:,sframe); type1_sframe_resp = type1_sframe_resp(:);
        type2_sframe_resp = -type2_resp(:,sframe); type2_sframe_resp = type2_sframe_resp(:);
        type3_sframe_resp = -type3_resp(:,sframe); type3_sframe_resp = type3_sframe_resp(:);
        type5_sframe_resp = -type5_resp(:,sframe); type5_sframe_resp = type5_sframe_resp(:);
        type6_sframe_resp = -type6_resp(:,sframe); type6_sframe_resp = type6_sframe_resp(:);
        type7_sframe_resp = -type7_resp(:,sframe); type7_sframe_resp = type7_sframe_resp(:);

        all_resp = ([type1_sframe_resp;type2_sframe_resp;type3_sframe_resp;type1_sframe_resp;type6_sframe_resp;type7_sframe_resp]-bias)/scaling;
        all_resp_cell = {(type3_sframe_resp-bias)/scaling;(type2_sframe_resp-bias)/scaling;(type1_sframe_resp-bias)/scaling;(type7_sframe_resp-bias)/scaling;(type6_sframe_resp-bias)/scaling;(type5_sframe_resp-bias)/scaling};

        subplot(1,3,2);
        % ax = axes;
        hold on;

        range = linspace(min(all_resp),max(all_resp),15);
        range = [2*range(1)-range(2),range,2*range(end)-range(end-1)];
        % range = min(all_resp):0.0015:max(all_resp);

        for i = 1:6
            EBdata = all_resp_cell{i};
            B = 5000;
            histfun = @(x) histcounts(x, range)/length(x);
            bootstrap_proportion = bootstrp(B, histfun, EBdata);
            %     sem_counts = std(bootstrap_counts);
            %     sem_proportions = sem_counts / length(EBdata);

            plotShadedErrorbar_flip_givenY((range(1:end-1)+range(2:end))/2,bootstrap_proportion,colors(i,:),0.3);
            yline(mean(EBdata,'omitnan'),'Color',colors(i,:),'Linestyle','--','Linewidth',1.5);

        end

        % legend({'','-0.12','','','-0.16','','','-0.2','','','0.2','','','0.16','','','0.2'});
        % legend boxoff
        set(gca,'FontSize',16,'LineWidth',0.5);
        % set(gcf,'Position',[100 500 800 400]);
        set(findall(gcf, 'Type', 'text'), 'FontName', 'Calibri'); % 所有文本
        % y_lim = ylim;
        % ylim([0 y_lim(2)]);
        % xlim(x_lim);
        % x_lim = [min(all_resp),max(all_resp)];
        ylim([prctile(all_resp,2) prctile(all_resp,98)]);
%         x_ticks = xticks;
        %
%         x_ticks_label = cellfun(@(x) num2str(x), num2cell(-x_ticks), 'UniformOutput', false);
%         xticklabels(x_ticks_label);
        
        subplot(1,3,3);
        hold on;
        i = region;
        types_mean= [mean(all_resp_cell{1}),mean(all_resp_cell{2}),mean(all_resp_cell{3}),mean(all_resp_cell{4}),mean(all_resp_cell{5}),mean(all_resp_cell{6})];

        plot(1:6,types_mean,'-','LineWidth',1.5,'Color',[0.7 0.7 0.7]);

        errorbar(1,mean(all_resp_cell{1}),sem(all_resp_cell{1}),'o',"MarkerSize",5,'LineWidth',1.5,'CapSize',10,'Color',colors(1,:),'MarkerFaceColor','auto');
        errorbar(2,mean(all_resp_cell{2}),sem(all_resp_cell{2}),'o',"MarkerSize",5,'LineWidth',1.5,'CapSize',10,'Color',colors(2,:),'MarkerFaceColor','auto');
        errorbar(3,mean(all_resp_cell{3}),sem(all_resp_cell{3}),'o',"MarkerSize",5,'LineWidth',1.5,'CapSize',10,'Color',colors(3,:),'MarkerFaceColor','auto');
        errorbar(4,mean(all_resp_cell{4}),sem(all_resp_cell{4}),'o',"MarkerSize",5,'LineWidth',1.5,'CapSize',10,'Color',colors(4,:),'MarkerFaceColor','auto');
        errorbar(5,mean(all_resp_cell{5}),sem(all_resp_cell{5}),'o',"MarkerSize",5,'LineWidth',1.5,'CapSize',10,'Color',colors(5,:),'MarkerFaceColor','auto');
        errorbar(6,mean(all_resp_cell{6}),sem(all_resp_cell{6}),'o',"MarkerSize",5,'LineWidth',1.5,'CapSize',10,'Color',colors(6,:),'MarkerFaceColor','auto');

        %        types_sem = [type3_sem(i,sframe),type2_sem(i,sframe),type1_sem(i,sframe),type7_sem(i,sframe),type6_sem(i,sframe),type5_sem(i,sframe)];
        box off;
        title(region_label(i+1).abbrev);
        set(gca,'Xcolor','none','FontSize',16);
        xlim([0 length(obs_type)]);
        % xticks = 0:length(obs_type);
        xticklabels({});
        %       xticklabels = {'{}','-0.12(L)','-0.16(L)','-0.2(L)','0.2(R)','0.16(R)','0.12(R)','{}'};
        % xticklabels = {'{}','L1','L2','L3','R3','R2','R1','{}'};
%         set(gca, 'XColor', 'none');
        % ylim(x_lim);
        sem_average = mean([sem(all_resp_cell{1}),sem(all_resp_cell{2}),sem(all_resp_cell{3}),sem(all_resp_cell{4}),sem(all_resp_cell{5}),sem(all_resp_cell{6})]);
        set(gcf,'Position',[100 500 880 230]);
        set(findall(gcf, 'Type', 'text'), 'FontName', 'Calibri'); % 所有文本
        ylim([-8*sem_average 8*sem_average]);
    end
end