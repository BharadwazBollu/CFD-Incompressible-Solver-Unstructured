close all
clear 
clc
%%

data = dlmread('u_centerline.dat');
y = data(:,2) ;
U = data(:,3) ;
data = dlmread('v_centerline.dat');
x = data(:,1) ;
V = data(:,4) ;

ghiaU=dlmread('ghiaU.dat');
choiU=dlmread('choiU.dat');
ghiaV=dlmread('ghiaV.dat');
choiV=dlmread('choiV.dat');

figure(1)
plot(U, y ,'r*',ghiaU(:,1),ghiaU(:,2),'kd',choiU(:,1),choiU(:,2),'ko')
xlabel(' U Velocity ');
ylabel(' Y - distance ');
title('U Velocity Center line')
legend('UPWIND','GHIA','CHOI','Location','northwest')
% print(gcf,'U_UPWIND_IMPLICIT.jpg','-dpng','-r300');
figure(2)
plot(x,V,'r*',ghiaV(:,1),ghiaV(:,2),'kd',choiV(:,1),choiV(:,2),'ko')
xlim([ 0 1 ] );
xlabel(' X - distance ');
ylabel(' V Velocity ');
title('V Velocity Center line ')
legend('UPWIND','GHIA','CHOI','Location','northeast')
% print(gcf,'V_UPWIND_IMPLICIT.jpg','-dpng','-r300');