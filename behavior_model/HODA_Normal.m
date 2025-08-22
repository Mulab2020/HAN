%% Map generation
x = 1:1000;
y = 1:1000;
[X,Y] = meshgrid(x,y);

%%
tic;
out_times = [];
for sig = [20 100 200]
out_times = [];
for k = 1:100
%     sig = 13;
    sigma=[sig^2 0;
           0 sig^2];
    rng(k);
    start_pos = [1000 1000;1 1;1 1000;1000 1;1 500; 500 1;1000 500; 500 1000];
    obs_num = 100;
    obstacle_pos = zeros(obs_num*size(start_pos,1),2);
    for j = 1:size(start_pos,1)
        obstacle_pos((j-1)*obs_num+1,:) = start_pos(j,:);
        for i = 2:obs_num
            if i > 2
                mu_diff = obstacle_pos((j-1)*obs_num+i-1,:)-obstacle_pos((j-1)*obs_num+i-2,:);
                P = mvnpdf([X(:),Y(:)],flip(obstacle_pos((j-1)*obs_num+i-1,:)+15*mu_diff/norm(mu_diff)),sigma);
            else
                P = mvnpdf([X(:),Y(:)],obstacle_pos((j-1)*obs_num+i-1,:),sigma);
            end
            P(sub2ind([1000 1000],obstacle_pos((j-1)*obs_num+i-1,1),obstacle_pos((j-1)*obs_num+i-1,2)))= 0;
            P_cum = cumsum(P);
            rand_pos = rand(1).*sum(P_cum(end),'all');
            [obstacle_pos((j-1)*obs_num+i,1),obstacle_pos((j-1)*obs_num+i,2)] = ind2sub([1000 1000],find(P_cum>rand_pos,1));
        end
    end
    params.step_size = 1;
    params.attend_dist = 11;
    params.attend_angle = 180;
    params.turning_angle = 5;
    params.noise = 5;
    params.attend_barrier_dist = 2;
    params.decay_tau = 30;
    init_heading_pool = [0:359];
    init_headings = zeros(2,length(init_heading_pool));
    for i = 1:length(init_heading_pool)
        init_headings(:,i) = [sin(init_heading_pool(i)/180*pi);cos(init_heading_pool(i)/180*pi)];
    end
    trajectories_no_bias = cell(1,length(init_heading_pool));
    trajectories_with_bias = cell(1,length(init_heading_pool));

    out_time = zeros(length(init_heading_pool),2);

    parfor i = 1:length(init_heading_pool)
        init_heading = init_headings(:,i);
        init_pos = [500 500];
        [trajectory_no_bias,~] = no_bias_fish_2(init_heading,init_pos,obstacle_pos,i,params,1000,1000,0);
        trajectories_no_bias{i} = trajectory_no_bias;
        [trajectory_with_bias,barrier_facing_pos,~] = with_bias_fish_2(init_heading,init_pos,obstacle_pos,i,params,1000,1000,0);
        trajectories_with_bias{i} = trajectory_with_bias;
        out_time(i,:) = [length(trajectory_no_bias) length(trajectory_with_bias)];
    end
    out_times = [out_times;out_time];
end
out_times_diffsigma{sig} = out_times;
end
toc;
%%
figure;
data = out_times_diffsigma{100};
bar(1,mean(data(:,1)),'FaceColor',[0.2, 0.2, 0.2],'LineWidth',2);
hold on;
bar(2,mean(data(:,2)),'FaceColor',[0.5, 0.5, 0.5],'LineWidth',2);

errorbar([1 2],mean(data),[0 0],[sem(data(:,1)) sem(data(:,2))],'Color','k','LineWidth',1,'CapSize',30,"LineStyle","none");
xticks([1 2]);
xticklabels({'no bias','with bias'});
box off
ylabel('Swim out steps');
set(gca,'Fontsize',16,'LineWidth',1);
%% Exploring efficiency
outMean_no_bias = zeros(1,100);
outMean_with_bias = zeros(1,100);
x_range = zeros(1,length(out_times_diffsigma));
for i =1:length(x_range)
    x_range(i) = isempty(out_times_diffsigma{i});
end
x_range = find(~x_range);
for i = x_range
    outMean_no_bias(i) = mean(out_times_diffsigma{i}(:,1));
    outMean_with_bias(i) = mean(out_times_diffsigma{i}(:,2));
end

figure;
outMean_no_bias(outMean_no_bias==0) = [];
outMean_with_bias(outMean_with_bias==0) = [];
outMean_no_bias = outMean_no_bias.^(-1);
outMean_with_bias = outMean_with_bias.^(-1);

plot(x_range,outMean_no_bias,'b','LineWidth',1);
hold on;
sem_1 = [];
for i = x_range
    sem_1 = [sem_1 sem(out_times_diffsigma{i}(:,1).^(-1))];
end
curve1 = outMean_no_bias + sem_1;
curve2 = outMean_no_bias - sem_1;
x = x_range;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0 0 1],'FaceAlpha',0.3,'EdgeColor','none');
plot(x_range,outMean_with_bias,'r','LineWidth',1);
sem_2 = [];
for i = x_range
    sem_2 = [sem_2 sem(out_times_diffsigma{i}(:,2).^(-1))];
end
curve1 = outMean_with_bias + sem_1;
curve2 = outMean_with_bias - sem_1;
x = x_range;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [1 0 0],'FaceAlpha',0.3,'EdgeColor','none');

% for i = x_range
%     if kstest2(out_times_diffsigma{i}(:,1),out_times_diffsigma{i}(:,2))
%         scatter(i,3.8e-4,15,'*','k');
%     end
% end
% plot([30, 50], [3.8e-4, 3.8e-4],'k', 'LineWidth',5);  
legend({'no bias','','with bias'}); legend boxoff;
set(gca,'FontSize',16,'LineWidth',1.2);
ylabel('Exploring efficiency (1/outTime)');
xlabel('\sigma');
% xlim([30 100]);
% ylim([3.1e-4 3.85e-4]);
box off;