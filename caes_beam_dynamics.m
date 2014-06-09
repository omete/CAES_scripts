% CAES Beam Dynamics
%% Variables:
%  time stdx stdy stdz avgz avgG avgr stdG nemixrms nemiyrms
% Command line import 
close all;
clear;
clc;
numfile = 5;
mode = 1;    % 1 for beam dynamics 2 for (Q and numpar)

% Beam dynamics
for i=1:numfile
fileID = fopen(['avgs' num2str(i) '.txt']);
if (mode == 1)
    formatSpec = '%f %f %f %f %f %f %f %f %f %f';
end
if (mode == 2)
    formatSpec = '%f %f %f';
end
DATA(i).Data = textscan(fileID, formatSpec,'HeaderLines',1);
fclose(fileID);
end

% Assign variables
if (mode == 1),
for i=1:numfile
%DATA(i).time = DATA(i).Data{1};
DATA(i).pos = DATA(i).Data{1};
DATA(i).stdx = DATA(i).Data{2};
DATA(i).stdy = DATA(i).Data{3};
DATA(i).stdz = DATA(i).Data{4};
DATA(i).avgZ = DATA(i).Data{5};
DATA(i).avgG = DATA(i).Data{6};
DATA(i).avgr = DATA(i).Data{7};
DATA(i).stdG = DATA(i).Data{8};
DATA(i).nemixrms = DATA(i).Data{9};
DATA(i).nemiyrms = DATA(i).Data{10};
end
end

if (mode == 2),
for i=1:numfile
DATA(i).time = DATA(i).Data{1};
%DATA(i).pos = DATA(i).Data{1};
DATA(i).Q = DATA(i).Data{2};
DATA(i).numpar = DATA(i).Data{3};
end 
end

%% Emittance at a particular location

for i=1:numfile
    emitt_poi(i) = DATA(i).nemixrms(find(DATA(i).pos == 6.5e-2));
end

emitt_poi_mean = round( mean(emitt_poi)*1e10)/1e10;
emitt_poi_std  = round( std(emitt_poi)*1e10)/1e10;


% Figures
figs = 0;  % 0 - do not save 1 - save
colours = jet(numfile);

%figure('units','normalized','outerposition',[0 0 1 1]) % for fullscreen
figure(1)
plot(DATA(1).pos(1)*1e3, DATA(1).nemixrms(1)*1e6,'.')
hold on;
for i=1:numfile
   h(i) = plot(DATA(i).pos*1e3,DATA(i).nemixrms*1e6,'-o','color',colours(i,:));
end
%Scheme 1
%stem([-10 10 13 150], max(DATA(1).nemixrms/2)*ones(1,4)*1e6,'m','linewidth',1) 
%Scheme 2
stem([-15 15 35 65], 20*ones(1,4),'m','linewidth',2) 
hold off;
xlabel('Position (mm)','fontsize',14)
ylabel('\epsilon (mm mrad)','fontsize',14)
xlim([-20 70])
ylim([0 12])
legend([h(1) h(2) h(3) h(4) h(5)], '1V','10V','100V','1kV','10kV');
%title(['<\epsilon> = ' num2str(emitt_poi_mean) ' m-rad, \sigma{\epsilon} = ' num2str(emitt_poi_std) ' m-rad']);
if (figs == 1)
saveas(gca, 'emittance.eps','epsc');
end


%figure('units','normalized','outerposition',[0 0 1 1])
figure(2)
h1 = plot(DATA(1).pos*1e3, DATA(1).stdx*1e3,'-b')
hold on;
for i=1:numfile
    h(i)=plot(DATA(i).pos*1e3, DATA(i).stdx*1e3,'-o','color',colours(i,:));
end
%Scheme1
%stem([-10 10 13 150], max(DATA(1).stdx)*ones(1,4)*1e3,'m','linewidth',1)
%Scheme2
stem([-15 15 35 65], max(DATA(1).stdx)*ones(1,4)*1e3,'m','linewidth',1)
xlabel('Position (mm)','fontsize',14)
ylabel('\sigma_{x,y} (mm)','fontsize',14)
xlim([-20 65])
%h=legend([h1 h2 h3 h4],'\sigma_x mesh','\sigma_x hole', '<r> mesh','<r> hole');
%leg = findobj(h,'type','text');
%set(leg, 'fontsize',14);
%legend([h(1) h(2) h(3)],'50ns','200ns','1\mus')
if (figs == 1)
saveas(gca, 'sigmaxr.eps','epsc');
end

t_inj = [100 200 300 400 500 600 700 800 900 1000 1200 1400 1600 1800 2000]

% Model: Gaussian distribution 
F = @(x,xdata)x(1)*xdata.^-2+x(2)*xdata+x(3);
%F = @(x,xdata)x(1)./exp(x(2).*xdata.^2) + x(3);
% Initial values
% Perform the fittings  
x0 =[1 0 1];
[x1,resnorm1,~,exitflag1,output1] = lsqcurvefit(F,x0,t_inj(3:15),emitt_poi(3:15));


figure(3)
semilogy(t_inj(3:15),emitt_poi(3:15), '-o')
hold on;
errorbar(t_inj, ones(length(t_inj),1)*2.4e-9,ones(length(t_inj),1)*0.2e-9)
semilogy(t_inj(),F(x1,t_inj()),'r')
hold off;
xlabel('Injection Time (ns)','fontsize',14)
ylabel('\epsilon (mm mrad)','fontsize',14)
xlim([0 2100])
%%
figure('units','normalized','outerposition',[0 0 1 1])
h1 = plot(DATA(1).pos*1e3, (DATA(1).stdz/3e8)*1e12)
hold on;
for i=1:numfile
    h(i)=plot(DATA(i).pos*1e3, (DATA(i).stdz/3e8)*1e12,'-o','color',colours(i,:))
end
stem([-10 10 13 150], max(DATA(1).stdz/3e8)*ones(1,4)*1e12/2,'m','linewidth',1)
hold off;
xlabel('s, position (mm)','fontsize',14)
ylabel('\sigma_{z} (ps)','fontsize',14)
%h=legend([h1 h2],'\sigma_z mesh','\sigma_z hole')
%leg = findobj(h,'type','text');
%set(leg, 'fontsize',14);
legend([h(1) h(2) h(3) h(4) h(5)],'1V','10V','100V','1kV','10kV')
xlim([-15 200])
if (figs == 1)
saveas(gca, 'sigmaz.eps','epsc');
end
%%
figure(4)
h1 = errorbar(DATA(1).pos(1)*1e3, DATA(1).avgG(1),DATA(1).stdG(1)/2)
hold on;
for i=1:numfile
    errorbar(DATA(i).pos*1e3, DATA(i).avgG,DATA(i).stdG/2,'-o','color',colours(i,:))
end
stem([-10 10 13 150], max(DATA(1).avgG)*ones(1,4),'m','linewidth',1)
hold off;
xlabel('s, position (mm)','fontsize',14)
ylabel('<\gamma>','fontsize',14)
%xlim([-15 500])
%ylim([1 1.025])
%h=legend([h1 h2],'\gamma mesh','\gamma hole');
%leg = findobj(h,'type','text');
%set(leg, 'fontsize',14);
colorbar;
if (figs == 1)
saveas(gca, 'gamma.eps','epsc');
end

figure(5)
h1 = plot(DATA(1).pos*1e3, DATA(1).stdG./DATA(1).avgG(1))
hold on;
for i=1:numfile
    plot(DATA(i).pos*1e3, DATA(i).stdG./DATA(i).avgG,'-o','color',colours(i,:))
end
stem([-10 10 13 150], max(DATA(1).stdG./DATA(1).avgG)*ones(1,4),'m','linewidth',1)
hold off;
xlabel('s, position (mm)','fontsize',14)
ylabel('\Delta\gamma / \gamma','fontsize',18)
%h=legend([h1 h2],'\sigma_z mesh','\sigma_z hole')
%leg = findobj(h,'type','text');
%set(leg, 'fontsize',14);
colorbar;
%xlim([-15 500])
if (figs == 1)
saveas(gca, 'delta_gamma.eps','epsc');
end

%% Plots at a particular locations

delta_volt = [1 2 3 4 5 6 7 8 9 10 11 12 13 14]; % Voltage difference between plate 1 and 2

%indx = find(DATA(1).pos == 0.1503);
indx = 221;
% Assign variables
for i=1:length(delta_volt)
pos_(i)      = DATA(i).pos(indx);
stdx_(i)     = DATA(i).stdx(indx);
stdy_(i)     = DATA(i).stdy(indx);
stdz_(i)     = DATA(i).stdz(indx);
avgZ_(i)     = DATA(i).avgZ(indx);
avgG_(i)     = DATA(i).avgG(indx);
avgr_(i)     = DATA(i).avgr(indx);
stdG_(i)     = DATA(i).stdG(indx);
nemixrms_(i) = DATA(i).nemixrms(indx);
nemiyrms_(i) = DATA(i).nemiyrms(indx);
end

deltaGG = stdG_./avgG_;

figure(1)
[hax, hline1, hline2] = plotyy(delta_volt, nemixrms_*1e6, delta_volt, deltaGG)
xlabel('\DeltaV (kV)','fontsize',18)
ylabel(hax(1),'\epsilon_{n,x} (mm mrad)','fontsize',18)
ylabel(hax(2),'\Delta\gamma / \gamma','fontsize',18)
set(hline1,'linestyle','o')
set(hline2,'linestyle','o')


%% Q and numpar

figure(1)
plot(DATA(1).time, DATA(1).Q*1e12,'r')
hold on;
plot(DATA(2).time, DATA(2).Q*1e12,'b')
plot(DATA(3).time, DATA(3).Q*1e12,'g')
hold off;
xlabel('Position (m)','fontsize',14)
ylabel('Q (pC)','fontsize',14)
legend('50ns','200ns','1\mus')

figure(2)
plot(DATA(1).time, DATA(1).numpar,'r')
hold on;
plot(DATA(2).time, DATA(2).numpar,'b')
plot(DATA(3).time, DATA(3).numpar,'g')
hold off;






