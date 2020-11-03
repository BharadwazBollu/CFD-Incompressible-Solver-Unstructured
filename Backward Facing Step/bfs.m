close all
clear
clc
%%
gartling7   = dlmread('gartling_bfs_X7.dat')    ;
gartling15  = dlmread('gartling_bfs_X15.dat')   ;

bfs7        = dlmread('bfs_X7.dat')    ;
bfs15       = dlmread('bfs_X15.dat')   ;

bfs7(:,1)   = bfs7(:,1)  - 0.5  ;
bfs15(:,1)  = bfs15(:,1) - 0.5  ;

hold on
plot(gartling7(:,1),gartling7(:,2),'r*','MarkerSize',12,'LineWidth',2);
plot(bfs7(:,2),bfs7(:,1),'k','LineWidth',2);
plot(gartling15(:,1),gartling15(:,2),'b*','MarkerSize',12,'LineWidth',2);
plot(bfs15(:,2),bfs15(:,1),'k--','LineWidth',2);
xlabel(' U ','FontSize',18,'FontWeight','bold');
ylabel(' Y ','FontSize',18,'FontWeight','bold');
set(gca,'FontSize',18,'FontWeight','bold');
legend('gartling x = 7 ','bfs x = 7 ','gartling x = 15 ','bfs x = 15 ','Location','northeast','FontSize',20,'FontWeight','bold');
axis([-0.1 1.2 -0.55 0.55])
% print(gcf,'bfs.jpg','-dpng','-r300');

