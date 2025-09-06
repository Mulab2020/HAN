function [trajectory,barrier_facing_pos] = no_bias_fish_2(init_heading,init_pos,obstacle_pos,rand_seed,params,x_lim,y_lim,plot_or_not)
curr_heading = init_heading;
curr_pos = init_pos;
trajectory = curr_pos;
barrier_facing_pos = [];
collide_barrier = false;
rng(rand_seed);
while true
    xy_diff = obstacle_pos-curr_pos;
    obs_dist = sqrt(xy_diff(:,1).^2 + xy_diff(:,2).^2);
    front_obs = (xy_diff./sqrt(xy_diff(:,1).^2+xy_diff(:,2).^2)*curr_heading)>=cos(params.attend_angle/2/180*pi);
    front_obs_dist = obs_dist(front_obs,:);

    L_obs_idx = obstacle_pos(:,2)*curr_heading(1)-curr_heading(2)*obstacle_pos(:,1)+curr_heading(2)*curr_pos(1)-curr_pos(2)*curr_heading(1)>0;
    R_obs_idx = obstacle_pos(:,2)*curr_heading(1)-curr_heading(2)*obstacle_pos(:,1)+curr_heading(2)*curr_pos(1)-curr_pos(2)*curr_heading(1)<0;
    L_obs = obstacle_pos(L_obs_idx,:);
    R_obs = obstacle_pos(R_obs_idx,:);
    if ~isempty(L_obs)&&~isempty(R_obs)
        [~,L_closest_idx] = min(obs_dist(L_obs_idx));
        L_closest_obs = L_obs(L_closest_idx,:);
        [~,R_closest_idx] = min(obs_dist(R_obs_idx));
        R_closest_obs = R_obs(R_closest_idx,:);
        barrier_direction = R_closest_obs-L_closest_obs;% 得到从左指向右的障碍物墙方向
        barrier_direction = barrier_direction/norm(barrier_direction);
        % 判断障碍物在左还是在右 y+ax+b>0 为左
        s1 = sqrt((R_closest_obs(1)-L_closest_obs(1))^2+(R_closest_obs(2)-L_closest_obs(2))^2);
        s2 = sqrt((L_closest_obs(1)-curr_pos(1))^2+(L_closest_obs(2)-curr_pos(2))^2);
        s3 = sqrt((R_closest_obs(1)-curr_pos(1))^2+(R_closest_obs(2)-curr_pos(2))^2);
        p = (s1+s2+s3)/2; % 计算三角形周长的一半
        barrier_dist = sqrt(p*(p-s1)*(p-s2)*(p-s3))*2/s1; % 计算现在位置到障碍墙之间的距离
        curr_2_obs = curr_pos-L_closest_obs;
        if barrier_dist < params.attend_barrier_dist && barrier_direction(1)*curr_2_obs(2)-barrier_direction(2)*curr_2_obs(1)<0  && s1<25% && abs((barrier_direction*curr_heading)/norm(barrier_direction)/norm(curr_heading))<0.5 % 两个障碍物在朝向的两侧，障碍物墙离动物足够近
            if barrier_direction*curr_heading > 0 % 障碍物墙应该向右躲避
                curr_heading = barrier_direction';
            else
                curr_heading = -barrier_direction';
            end
            barrier_facing_pos = [barrier_facing_pos;curr_pos];
            collide_barrier = true;
        end
    end
    if min(front_obs_dist)<params.attend_dist & ~collide_barrier
        [~,closest_obs_idx] = min(front_obs_dist);
        front_obs_pos = obstacle_pos(front_obs,:);
        closest_obs = front_obs_pos(closest_obs_idx,:);
        closest_obs_LR = closest_obs(2)*curr_heading(1)-curr_heading(2)*closest_obs(1)+curr_heading(2)*curr_pos(1)-curr_pos(2)*curr_heading(1);
        obs_theta = acos((closest_obs-curr_pos)*curr_heading/norm(closest_obs-curr_pos)/norm(curr_heading));
        scaling_factor = 4*exp(-sum((closest_obs-curr_pos).^2)/200); % *exp(-obs_theta/params.scaling_factor_tau1)
        if closest_obs_LR > 0 % && length(trajectory)-last_barrier_turn_time>2 || (length(trajectory)-last_barrier_turn_time<=2 && last_barrier_turn > 0) % 障碍物在左侧
            curr_heading = [cos(params.turning_angle/180*pi*scaling_factor) sin(params.turning_angle/180*pi*scaling_factor);
                -sin(params.turning_angle/180*pi*scaling_factor) cos(params.turning_angle/180*pi*scaling_factor)]*curr_heading;
        else
            curr_heading = [cos(params.turning_angle/180*pi*scaling_factor) -sin(params.turning_angle/180*pi*scaling_factor);
                sin(params.turning_angle/180*pi*scaling_factor) cos(params.turning_angle/180*pi*scaling_factor)]*curr_heading;
        end
    end
    if ~collide_barrier
        rand_degree_rad = params.noise*randn/180*pi; % 添加每次游泳的随机噪音
        noise = [cos(rand_degree_rad) sin(rand_degree_rad);
            -sin(rand_degree_rad) cos(rand_degree_rad)];
        curr_heading = noise*curr_heading;
    else
        collide_barrier = false;
    end
    last_pos = curr_pos;
    curr_pos = curr_pos+curr_heading'*params.step_size;
    trajectory = [trajectory;curr_pos];
%     if length(trajectory) == 3120
%         continue;
%     end
    if (curr_pos(1)<0 || curr_pos(1)>x_lim || curr_pos(2)<0 || curr_pos(2)>y_lim) || length(trajectory)==50000
        break;
    end
end
if plot_or_not
    figure;
    hold on;
    scatter(obstacle_pos(:,1),obstacle_pos(:,2),20,'filled');
    axis equal;
    plot(trajectory(:,1),trajectory(:,2));
    box on
    xlim([0 x_lim]);
    ylim([0 y_lim]);
    set(gcf,'position',[1241   777   807   749]);
    disp(['no bias trajectory length: ',num2str(length(trajectory))]);
end
end