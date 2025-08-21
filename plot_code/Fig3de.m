load('Fig3de_data_CD_intvl_resp.mat');
%%
set(0, 'DefaultFigureColor', 'w');
Rstim_range = 7:27;
Lstim_range = 6:26;
frame_rate = (length(Rstim_range)-1)/8;
L_control_color = 'r';
R_control_color = 'b';
L_aftPerturb_color = 'm';
R_aftPerturb_color = 'c';
figure;

subplot(2,1,1);
hold on;

plotShadedErrorbar_dashed(L_control_intvl_resp_all(:,Lstim_range),frame_rate,L_control_color,0.2);
plotShadedErrorbar_dashed(R_control_intvl_resp_all(:,Lstim_range),frame_rate,R_control_color,0.2);
plotShadedErrorbar(L_perturb_intvl_resp_all(L_stay_trials,Lstim_range),frame_rate,L_aftPerturb_color,0.2);
plotShadedErrorbar(L_perturb_intvl_resp_all(L_switch_trials,Lstim_range),frame_rate,R_aftPerturb_color,0.2);

ylabel('CD projection (a.u.)');
xlabel('Time (s)');

legend({'','ctrl L','','ctrl R','','stay L','','switch L'});
legend boxoff
set(gca,'FontSize',16,'LineWidth',1.2);

subplot(2,1,2);
hold on;

plotShadedErrorbar_dashed(L_control_intvl_resp_all(:,Rstim_range),frame_rate,L_control_color,0.2);
plotShadedErrorbar_dashed(R_control_intvl_resp_all(:,Rstim_range),frame_rate,R_control_color,0.2);
plotShadedErrorbar(R_perturb_intvl_resp_all(R_stay_trials,Rstim_range),frame_rate,R_aftPerturb_color,0.2);
plotShadedErrorbar(R_perturb_intvl_resp_all(R_switch_trials,Rstim_range),frame_rate,L_aftPerturb_color,0.2);

ylabel('CD projection (a.u.)');
xlabel('Time (s)');

legend({'','ctrl L','','ctrl R','','stay R','','switch R'});
legend boxoff
set(gca,'FontSize',16,'LineWidth',1.2);

set(gcf,'Position',[1300 800 800 1000]);
%%

R_switch_intvl_resp_all = R_perturb_intvl_resp_all(R_switch_trials,:);
L_switch_intvl_resp_all = L_perturb_intvl_resp_all(L_switch_trials,:);
R_stay_intvl_resp_all = R_perturb_intvl_resp_all(R_stay_trials,:);
L_stay_intvl_resp_all = L_perturb_intvl_resp_all(L_stay_trials,:);

L_switch_resp = L_switch_intvl_resp_all(:,Lstim_range);
L_stay_resp = L_stay_intvl_resp_all(:,Lstim_range);
R_switch_resp = R_switch_intvl_resp_all(:,Rstim_range);
R_stay_resp = R_stay_intvl_resp_all(:,Rstim_range);
L_ctrl_resp = L_control_intvl_resp_all(:,Lstim_range);
R_ctrl_resp = R_control_intvl_resp_all(:,Lstim_range);

figure;
L_control_color = 'r';
R_control_color = 'b';
L_aftPerturb_color = 'm';
R_aftPerturb_color = 'c';

final_range = (24:26)-5;
end_points = [mean(R_ctrl_resp(:,final_range),'all','omitnan') mean(L_ctrl_resp(:,final_range),'all','omitnan')];
bias = mean(end_points);

% curr_ctrl_resp = [R_ctrl_resp-bias;-(L_ctrl_resp-bias)];
% contra_ctrl_resp = -curr_ctrl_resp;
% switch_resp = [R_switch_resp-bias;-(L_switch_resp-bias)];
% stay_resp = [R_stay_resp-bias;-(L_stay_resp-bias)];

curr_ctrl_resp = [R_ctrl_resp;-L_ctrl_resp];
contra_ctrl_resp = -curr_ctrl_resp;
switch_resp = [R_switch_resp;-L_switch_resp];
stay_resp = [R_stay_resp;-L_stay_resp];


hold on;
plotShadedErrorbar_dashed(curr_ctrl_resp,frame_rate,L_control_color,0.2);
plotShadedErrorbar_dashed(contra_ctrl_resp,frame_rate,R_control_color,0.2);
plotShadedErrorbar(switch_resp,frame_rate,R_aftPerturb_color,0.2);
plotShadedErrorbar(stay_resp,frame_rate,L_aftPerturb_color,0.2);
% xline(1/frame_rate);
set(gca,'FontSize',16,'LineWidth',1.2);
legend({'','current state','','contra state','','switch','','stay'});
legend boxoff
ylabel('CD projection (a.u.)');
xlabel('Time (s)')
set(gcf,'Position',[1300 800 800 500]);
%%
bin_size = 0.27;
L_control_color = 'r';
R_control_color = 'b';
L_aftPerturb_color = 'm';
R_aftPerturb_color = 'c';
final_range = 17:21;
% curr_control_state = (curr_ctrl_resp(:,final_range)+curr_ctrl_resp(:,final_range+1)+curr_ctrl_resp(:,final_range-1))/3;
% curr_control_state = curr_control_state(:);
% 
% contra_control_state = (contra_ctrl_resp(:,final_range)+contra_ctrl_resp(:,final_range+1)+contra_ctrl_resp(:,final_range-1))/3;
% contra_control_state = contra_control_state(:);
% switch_state = (switch_resp(:,final_range)+switch_resp(:,final_range-1)+switch_resp(:,final_range+1))/3;
% switch_state = switch_state(:);
% stay_state = (stay_resp(:,final_range)+stay_resp(:,final_range-1)+stay_resp(:,final_range+1))/3;
% stay_state = stay_state(:);
curr_control_state = curr_ctrl_resp(:,final_range);
curr_control_state = curr_control_state(:);

contra_control_state = contra_ctrl_resp(:,final_range);
contra_control_state = contra_control_state(:);
switch_state = switch_resp(:,final_range);
switch_state = switch_state(:);
stay_state = stay_resp(:,final_range);
stay_state = stay_state(:);
figure;
ax = axes;
yyaxis(ax, 'left'); 
ax.YColor = 'k';
ylabel('Histogram probability');
hold on;
histogram(curr_control_state,'FaceColor',L_control_color,'BinWidth',bin_size,'Normalization','probability','EdgeColor','none','FaceAlpha',0.5);
histogram(contra_control_state,'FaceColor',R_control_color,'BinWidth',bin_size,'Normalization','probability','EdgeColor','none','FaceAlpha',0.5);
histogram(stay_state,'BinWidth',bin_size,'Normalization','probability','EdgeColor',R_aftPerturb_color,'DisplayStyle','stairs','LineWidth',1.5);
histogram(switch_state,'BinWidth',bin_size,'Normalization','probability','EdgeColor',L_aftPerturb_color,'DisplayStyle','stairs','LineWidth',1.5);
set(gca,'FontSize',16,'LineWidth',1.2);
yyaxis(ax, 'right'); 
all_states = [curr_control_state(:);contra_control_state(:);stay_state(:);switch_state(:)];
x = linspace(min(all_states), max(all_states), 100);
pd = fitdist(curr_control_state(:), 'Normal');
y = pdf(pd, x);
plot(x, y, L_control_color, 'LineWidth', 1.5); % 绘制拟合曲线
pd = fitdist(contra_control_state(:), 'Normal');
y = pdf(pd, x);
plot(x, y, R_control_color, 'LineWidth', 2,'LineStyle','-'); % 绘制拟合曲线
pd = fitdist(switch_state(:), 'Normal');
y = pdf(pd, x);
plot(x, y, L_aftPerturb_color, 'LineWidth', 2,'LineStyle','--'); % 绘制拟合曲线
pd = fitdist(stay_state(:), 'Normal');
y = pdf(pd, x);
plot(x, y, R_aftPerturb_color, 'LineWidth', 2,'LineStyle','--'); % 绘制拟合曲线
legend({'current state','contra state ','stay','switch'});
legend boxoff
ylabel('Fitted probability density');
ax.YColor = 'k';
set(gcf,'Position',[1300 800 800 500]);