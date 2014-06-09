%% PLOTS for the extraction region
%% Filed output from SUPERFISH
% Data List
%      X             Y              Ex            Ey            |E|           V 
%     (cm)          (cm)           (V/cm)        (V/cm)        (V/cm)         (V) 
data_field = load('field01.txt');
X  = data_field(:,1);
Y  = data_field(:,2);
EX = data_field(:,3);
EY = data_field(:,4);
aE = data_field(:,5);
V  = data_field(:,6);

figure(1) % Axial electric field EX
h2 = plot(X, V, '-or')
hold on;
h3 = plot([6.5 6.5],[-20000 5000],'--m')
plot([8.5 8.5],[-20000 5000],'--m')
plot([8.8 8.8],[-20000 5000],'--m')
h1 = plot(X, EX, '-ob')
hold off;
xlabel('Position (cm)','fontsize',14);
ylabel('EX (V/cm), V (V)','fontsize',14);
legend([h1 h2 h3], 'EX', 'V','Electrode locations')

