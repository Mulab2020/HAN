load('ExtDataFig6d_Bias_integrator_correlation.mat');
%%
figure;
bar_data = zeros(1,length(avoiding));
err_data = zeros(1,length(avoiding));
for i = 1:4
    bar_data(i) = mean(-avoiding{i});
    err_data(i) = sem(-avoiding{i});
end
errorbar(bar_data,err_data,'-o','MarkerSize',10,'CapSize',10,'LineWidth',2,'Color','k');
xlabel('Norm. {\Delta}F/F');
ylabel('Bias strength (a.u.)');
set(gcf,'Position',[100 100 500 500]);
set(gca,'Color','none','FontSize',16,'LineWidth',0.5, 'FontName', 'Calibri');
% set(gcf,'Position',[100 500 800 400]);
set(findall(gcf, 'Type', 'text'), 'FontName', 'Calibri'); % 所有文本
box off