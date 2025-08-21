function plotShadedErrorbar(data,frame_rate,color,shadealpha)
hold on;
sem_dev = std(data,1,'omitnan')./sqrt(sum(~isnan(data)));
y = mean(data,'omitnan');
curve1 = y + sem_dev;
curve2 = y - sem_dev;
nanIdx = isnan(curve1);

x = (0:size(data,2)-1)/frame_rate;
x = x(~nanIdx);
x2 = [x, fliplr(x)];

inBetween = [curve1(~nanIdx), fliplr(curve2(~nanIdx))];
fill(x2, inBetween, color,'FaceAlpha',shadealpha,'EdgeColor','none');
plot(x,y(~nanIdx),'LineWidth',1.5,'Color',color);
end








