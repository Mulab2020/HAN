function plotShadedErrorbar_flip_givenY(y,data,color,shadealpha)
hold on;
sem_dev = std(data,1,'omitnan');
x = mean(data,1,'omitnan');
curve1 = x + sem_dev;
curve2 = x - sem_dev;
nanIdx = isnan(curve1);


y2 = [y, fliplr(y)];
inBetween = [curve2(~nanIdx),fliplr(curve1(~nanIdx))];
fill(inBetween, y2, color,'FaceAlpha',shadealpha,'EdgeColor','none');
plot(x(~nanIdx),y,'LineWidth',2,'Color',color);
end