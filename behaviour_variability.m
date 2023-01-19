function behaviour_variability(session, Area,Nbins)
event=2;
Ndir=8;
colour_dir=hsv(Ndir);
ms=20;% this number doesn't matter because its related to the neural activity, not the behaviour
t_1=-0.5;
t_2=0.5;
from_to=[0.05 2]; %select movements between 50 ms and 2s
doplot=0;
colour_dur=plasma(Nbins);
[~,direction,~,~,~,~,dist_mov,mov_duration,max_speed]=neural_data_per_duration(session,Area,ms,t_1,t_2,event,from_to,doplot);

subplot(4,4,13)
h=polarhistogram(direction,Ndir);
Edges=h.BinEdges;
direction1=ceil(Ndir*(direction+pi)/(2*pi));
% make the bars coloured by direction
for i=1:Ndir
    polarhistogram(direction(direction1==i),Edges,'FaceColor',colour_dir(i,:));
    hold on
end
    
    
subplot(4,4,14)
hold on
h=histogram(mov_duration*1000,0:100:max(mov_duration)*1000,'normalization','probability','EdgeAlpha',0,'FaceColor','k');
middlebin=[250 350 450 550];
for i=1:numel(middlebin)
bar(middlebin(i),h.Values(i+2),100,'FaceColor',colour_dur(i,:))
end

box off
xlabel('Movement duration [ms]')
ylabel('Fraction of movements')
ylim([0 0.6])


subplot(4,4,15)
histogram(max_speed,0:1:max(max_speed),'normalization','probability','EdgeAlpha',0,'FaceColor','k')
box off
xlabel('Max speed [cm/s]')
ylim([0 0.1])

subplot(4,4,16)
histogram(dist_mov,0:1:max(dist_mov),'normalization','probability','EdgeAlpha',0,'FaceColor','k')
xlabel('Distance [cm]')
ylim([0 0.4])
box off

end