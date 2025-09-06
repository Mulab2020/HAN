load('ExtDataFig1b_turn_power.mat');
%%
CM = othercolor('PiYG4');
data1 = rmoutliers(All_R_trial_power);
data2 = rmoutliers(All_L_trial_power);
% data1 = All_last1_R_trial_power;
% data2 = All_last1_L_trial_power;


sorted_data = sort(data1);
n = length(data1);
ecdf_values = (1:n) / n;  % 经验 CDF的值

% 计算置信区间的宽度 (DKW 不等式)
alpha = 0.05;  % 95% 置信区间（alpha = 0.05）
epsilon = sqrt((log(2 / alpha)) / (2 * n));  % DKW不等式给出的置信带宽度

% 计算置信区间的上下限
upper_bound = ecdf_values + epsilon;  % 经验 CDF + 置信带宽度
upper_bound(upper_bound>1) = 1;
lower_bound = ecdf_values - epsilon;  % 经验 CDF - 置信带宽度
lower_bound(lower_bound<0) = 0;
% 绘制经验 CDF 和 置信区间
figure;
hold on;
yyaxis right;
plot(sorted_data, ecdf_values,'-', 'LineWidth', 2, 'Color', CM(end,:));  % 经验 CDF
% plot(sorted_data, upper_bound, '--r', 'LineWidth', 1.5);  % 置信区间上限
% plot(sorted_data, lower_bound, '--r', 'LineWidth', 1.5);  % 置信区间下限
fill([sorted_data,flip(sorted_data)],[lower_bound,flip(upper_bound)],CM(end,:),'FaceAlpha',0.2,'EdgeColor','none');

sorted_data = sort(data2);
n = length(data2);
ecdf_values = (1:n) / n;  % 经验 CDF的值
% 计算ECDF
[f1, x1] = ecdf(data1);
[f2, x2] = ecdf(data2);

% 去除重复x
[x1_unique, ia1] = unique(x1);
f1_unique = f1(ia1);

[x2_unique, ia2] = unique(x2);
f2_unique = f2(ia2);

% 合并x轴点并去重
x_common = unique([x1_unique; x2_unique]);

% 插值
f1_interp = interp1(x1_unique, f1_unique, x_common, 'previous', 'extrap');
f2_interp = interp1(x2_unique, f2_unique, x_common, 'previous', 'extrap');
diff_cdf = abs(f1_interp - f2_interp);
[max_diff, idx_max] = max(diff_cdf);
x_max = x_common(idx_max);
f1_max = f1_interp(idx_max);
f2_max = f2_interp(idx_max);
plot([x_max x_max], [f1_max, f2_max], 'k--', 'LineWidth', 1.5); % 连接线


% 计算置信区间的宽度 (DKW 不等式)
alpha = 0.05;  % 95% 置信区间（alpha = 0.05）
epsilon = sqrt((log(2 / alpha)) / (2 * n));  % DKW不等式给出的置信带宽度

% 计算置信区间的上下限
upper_bound = ecdf_values + epsilon;  % 经验 CDF + 置信带宽度
upper_bound(upper_bound>1) = 1;
lower_bound = ecdf_values - epsilon;  % 经验 CDF - 置信带宽度
lower_bound(lower_bound<0) = 0;
% 绘制经验 CDF 和 置信区间
plot(sorted_data, ecdf_values,'-', 'LineWidth', 2, 'Color', CM(1,:));  % 经验 CDF
% plot(sorted_data, upper_bound, '--b', 'LineWidth', 1.5);  % 置信区间上限
% plot(sorted_data, lower_bound, '--b', 'LineWidth', 1.5);  % 置信区间下限
fill([sorted_data,flip(sorted_data)],[lower_bound,flip(upper_bound)],CM(1,:),'FaceAlpha',0.2,'EdgeColor','none');
ylabel('Empirical CDF');
set(gca, 'YColor', 'k');

edges = linspace(min([data1,data2]), max([data1,data2]), 20);  % 设置 50 个区间
yyaxis left;
histogram(data1,edges, 'Normalization', 'probability','DisplayStyle','stairs','EdgeColor',CM(end,:),'LineWidth',1.5);
hold on;
histogram(data2,edges, 'Normalization', 'probability','DisplayStyle','stairs','EdgeColor',CM(1,:),'LineWidth',1.5);
% ylim([0 0.17]);

% 设置图形
% title('Empirical CDF with 95% Confidence Interval (DKW Inequality)');
xlabel('Turn power (a.u.)');
ylabel('Probability');
% legend('Empirical CDF', '95% Confidence Interval', 'Location', 'southeast');
% grid on;
legend({'R','L'});
legend boxoff;
set(gca,'FontSize',16,'LineWidth',1.2, 'YColor', 'k');
hold off;

%% Shuffle
CM = othercolor('PiYG4');
data1 = rmoutliers(All_R_trial_power);
data2 = rmoutliers(All_L_trial_power);
% data1 = All_last1_R_trial_power;
% data2 = All_last1_L_trial_power;
% Doing Shuffling
combined = [data1, data2];
shuffled_idx = randperm(length(combined));
shuffled_combined = combined(shuffled_idx);
data1 = shuffled_combined(1:length(data1));
data2 = shuffled_combined(1+length(data1):end);

sorted_data = sort(data1);
n = length(data1);
ecdf_values = (1:n) / n;  % 经验 CDF的值

% 计算置信区间的宽度 (DKW 不等式)
alpha = 0.05;  % 95% 置信区间（alpha = 0.05）
epsilon = sqrt((log(2 / alpha)) / (2 * n));  % DKW不等式给出的置信带宽度

% 计算置信区间的上下限
upper_bound = ecdf_values + epsilon;  % 经验 CDF + 置信带宽度
upper_bound(upper_bound>1) = 1;
lower_bound = ecdf_values - epsilon;  % 经验 CDF - 置信带宽度
lower_bound(lower_bound<0) = 0;
% 绘制经验 CDF 和 置信区间
figure;
hold on;
yyaxis right;
plot(sorted_data, ecdf_values,'-', 'LineWidth', 2, 'Color', CM(end,:));  % 经验 CDF
% plot(sorted_data, upper_bound, '--r', 'LineWidth', 1.5);  % 置信区间上限
% plot(sorted_data, lower_bound, '--r', 'LineWidth', 1.5);  % 置信区间下限
fill([sorted_data,flip(sorted_data)],[lower_bound,flip(upper_bound)],CM(end,:),'FaceAlpha',0.2,'EdgeColor','none');

sorted_data = sort(data2);
n = length(data2);
ecdf_values = (1:n) / n;  % 经验 CDF的值
% 计算ECDF
[f1, x1] = ecdf(data1);
[f2, x2] = ecdf(data2);

% 去除重复x
[x1_unique, ia1] = unique(x1);
f1_unique = f1(ia1);

[x2_unique, ia2] = unique(x2);
f2_unique = f2(ia2);

% 合并x轴点并去重
x_common = unique([x1_unique; x2_unique]);

% 插值
f1_interp = interp1(x1_unique, f1_unique, x_common, 'previous', 'extrap');
f2_interp = interp1(x2_unique, f2_unique, x_common, 'previous', 'extrap');
diff_cdf = abs(f1_interp - f2_interp);
[max_diff, idx_max] = max(diff_cdf);
x_max = x_common(idx_max);
f1_max = f1_interp(idx_max);
f2_max = f2_interp(idx_max);
plot([x_max x_max], [f1_max, f2_max], 'k--', 'LineWidth', 1.5); % 连接线


% 计算置信区间的宽度 (DKW 不等式)
alpha = 0.05;  % 95% 置信区间（alpha = 0.05）
epsilon = sqrt((log(2 / alpha)) / (2 * n));  % DKW不等式给出的置信带宽度

% 计算置信区间的上下限
upper_bound = ecdf_values + epsilon;  % 经验 CDF + 置信带宽度
upper_bound(upper_bound>1) = 1;
lower_bound = ecdf_values - epsilon;  % 经验 CDF - 置信带宽度
lower_bound(lower_bound<0) = 0;
% 绘制经验 CDF 和 置信区间
plot(sorted_data, ecdf_values,'-', 'LineWidth', 2, 'Color', CM(1,:));  % 经验 CDF
% plot(sorted_data, upper_bound, '--b', 'LineWidth', 1.5);  % 置信区间上限
% plot(sorted_data, lower_bound, '--b', 'LineWidth', 1.5);  % 置信区间下限
fill([sorted_data,flip(sorted_data)],[lower_bound,flip(upper_bound)],CM(1,:),'FaceAlpha',0.2,'EdgeColor','none');
ylabel('Empirical CDF');
set(gca, 'YColor', 'k');

edges = linspace(min([data1,data2]), max([data1,data2]), 20);  % 设置 50 个区间
yyaxis left;
histogram(data1,edges, 'Normalization', 'probability','DisplayStyle','stairs','EdgeColor',CM(end,:),'LineWidth',1.5);
hold on;
histogram(data2,edges, 'Normalization', 'probability','DisplayStyle','stairs','EdgeColor',CM(1,:),'LineWidth',1.5);
% ylim([0 0.17]);

% 设置图形
% title('Empirical CDF with 95% Confidence Interval (DKW Inequality)');
xlabel('Turn power (a.u.)');
ylabel('Probability');
% legend('Empirical CDF', '95% Confidence Interval', 'Location', 'southeast');
% grid on;
legend({'R','L'});
legend boxoff;
set(gca,'FontSize',16,'LineWidth',1.2, 'YColor', 'k');
hold off;