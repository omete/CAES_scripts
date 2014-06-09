%% CAES Beam Tracking Analyses
%% Variables:
%  x y rxy z G time
% Command line import 
close all;
clear;
clc;
numfile = 5;
for i=1:numfile
fileID = fopen(['traj' num2str(i) '.txt']);
formatSpec = '%f %f %f %f %f %f';
DATA(i).Data = textscan(fileID, formatSpec);
fclose(fileID);
end

% Assign variables
for i=1:numfile
DATA(i).x = DATA(i).Data{1};
DATA(i).y = DATA(i).Data{2};
DATA(i).rxy = DATA(i).Data{3};
DATA(i).z = DATA(i).Data{4};
DATA(i).G = DATA(i).Data{5};
DATA(i).time = DATA(i).Data{6};
end
%%
% Select a longitudinal position
file = 6;

% Correction
%indc = find(DATA(file).time == 9.995e-11);
%DATA(file).time(indc) = 0;

for i=1:350

%slice_ini = find(DATA(file).time == -4e-14);
slice_final = find(round(DATA(file).time*1e8*10)/(1e8*10) == (i)*1e-9);

figure(1)
plot3(DATA(file).z(slice_final),DATA(file).x(slice_final),DATA(file).y(slice_final),'ob')
%plot3(DATA(file).z(slice_ini),DATA(file).x(slice_ini),DATA(file).y(slice_ini),'ob')
hold on
%plot3(DATA(file).z(slice_final),DATA(file).x(slice_final),DATA(file).y(slice_final),'or')
%plot3(-10e-3,0,0,'o','color',[0.5 0.5 0.5],'linewidth',20)
%plot3(10e-3,0,0,'o','color',[0.5 0.5 0.6],'linewidth',20)
%plot3(13e-3,0,0,'o','color',[0.5 0.5 0.7],'linewidth',20)

plot3(-15e-3,0,0,'o','color',[0.5 0.5 0.6],'linewidth',20)
plot3(15e-3,0,0,'o','color',[0.5 0.5 0.7],'linewidth',20)
plot3(35e-3,0,0,'o','color',[0.5 0.5 0.7],'linewidth',20)
plot3(65e-3,0,0,'o','color',[0.5 0.5 0.8],'linewidth',20)
hold off
xlim([-20e-3 80e-3])
%xlim([-15e-3 50e-3]);
%xlim([-0.3e-6 0.3e-6]);
pbaspect([10 1 1])
%view(0,90)
xlabel('z')
ylabel('x')
zlabel('y')
ylim([-2e-3 2e-3]);
zlim([-2e-3 2e-3]);
grid on;
title([{'Scheme 2, Q = 1pC, Scheff on, E_{ini}=300meV, Electrodes (0, 0, +1V), CP (+200V)'}, {'Frame = ' num2str(i)}])
%title(['Scheme 1, Q = 1pC, Scheff on, Electrodes -15kV, -5kV, 0 --- Frame = ' num2str(i)])
%saveas(gca, ['xyz ' num2str(i) '.eps'],'epsc');
F(i) = getframe(gcf);
%pause(1)
end
%my_mov = movie(F);
movie2avi(F, ['Scheme2__scheff_sideview' num2str(file) '.avi'],'compression','none')
clear F;
%%
% REAL SPACE
figure(1)
plot(DATA(file).x(slice_ini),DATA(file).y(slice_ini),'o')

% HISTOS
numbins1 = 50;
figure(2);
subplot(1,2,1)
h1 = histfit(DATA(file).x(slice_ini),numbins1);
xlabel('x (m)')
ylabel('# of particles')
subplot(1,2,2)
h1 = histfit(DATA(file).y(slice_ini),numbins1);
xlabel('y (m)')
ylabel('# of particles')


%% Time evolution x-y-z
file = 1;
slice_ini = find(DATA(file).time == 0);
slice_final = find(DATA(file).time == 0*1e-3);
% 3D TRACKS
figure(3)
plot3(DATA(file).z(slice_ini),DATA(file).x(slice_ini),DATA(file).y(slice_ini),'ob')
hold on;
plot3(DATA(file).z(slice_final),DATA(file).x(slice_final),DATA(file).y(slice_final),'or')
hold off;
%xlim([-0.001 0.003])
%legend('initial','final')
pbaspect([10 1 1])
xlabel('z')
ylabel('x')
zlabel('y')
grid on;
saveas(gcf,['shot_' num2str(slice_final*1e7) '.eps'],'epsc')

%% Time evolution x-z plane
file = 6;

% Correction
%indc = find(DATA(file).time == 3.155e-30);
%DATA(file).time(indc) = 0;

target = 0*1e-9;
t = find(DATA(file).time == target);

% 2D TRACKS
figure(3)
plot(DATA(file).z(t),DATA(file).x(t),'ob')
xlim([0 300e-3]);
ylim([-2e-3 2e-3]);
xlabel('z')
ylabel('x')
title(['Time = ' num2str(target)]);
grid on;
%saveas(gcf,['shot_' num2str(file) '_' num2str(target) '.eps'],'epsc')

%% Time evolution x-y plane

file = 6;
t = find(DATA(file).time == 1*1e-6);

% 2D TRACKS
figure(3)
plot(DATA(file).x(t),DATA(file).y(t),'+b')
%xlim([-0.5e-6 0.5e-6])
%legend('initial','final')
xlabel('y','fontsize',14)
ylabel('x','fontsize',14)
grid on;
saveas(gcf,['profile_' num2str(file) '_' num2str(t(1)) '.eps'],'epsc')

%% Positions evolution x-y plane (BETTER)

file = 6;
p = find(DATA(file).time ==64*1e-3);

% 2D TRACKS
figure(3)
plot(DATA(file).x(p)*1e3,DATA(file).y(p)*1e3,'+b')
xlim([-4.5 4.5])
ylim([-4.5 4.5])
%legend('initial','final')
xlabel('y (mm)','fontsize',14)
ylabel('x (mm)','fontsize',14)
grid on;
saveas(gcf,['profile_' num2str(file) '_' num2str(p(1)) '.eps'],'epsc')

%% Draw all trajectories
file =6;
% 3D TRACKS

r = 52e-3;
rh = 10e-3;
theta = -pi:0.01:pi;
x = r*cos(theta);
y = r*sin(theta);
xh = rh*cos(theta);
yh = rh*sin(theta);
%figure(1)
%plot3(zeros(1,numel(x)),x,y,'b')


figure(1)
plot3(DATA(file).z(1),DATA(file).x(1),DATA(file).y(1),'+r')
hold on;
% Scheme 1
%fill3(ones(1,numel(x))*-10e-3,x,y,'b')
%fill3(ones(1,numel(x))*10e-3,x,y,'b')
%fill3(ones(1,numel(x))*10e-3,xh,yh,'w')
%fill3(ones(1,numel(x))*13e-3,x,y,'b')
%fill3(ones(1,numel(x))*65e-3,x,y,'g')
% Scheme 2
fill3(ones(1,numel(x))*-15e-3,x,y,'b')
fill3(ones(1,numel(x))*15e-3,x,y,'b')
fill3(ones(1,numel(x))*15e-3,xh,yh,'w')
fill3(ones(1,numel(x))*35e-3,x,y,'b')
fill3(ones(1,numel(x))*65e-3,x,y,'g')
for i=1:length(DATA(file).time)
plot3(DATA(file).z(i),DATA(file).x(i),DATA(file).y(i),'+r');
end
hold off;
xlim([-20e-3 80e-3]);   % z axis
ylim([-150e-3 60e-3]);   % x axis 
zlim([-60e-3 60e-3]);   % y axis
%pbaspect([3 1 1])
%view(3)
view(0,90)
xlabel('z (m)','fontsize',14)
ylabel('x (m)','fontsize',14)
zlabel('y (m)','fontsize',14)
grid on;









