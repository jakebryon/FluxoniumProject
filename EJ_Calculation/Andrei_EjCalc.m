%% estimate Ej from RT resistance vs area
% clc;
charge = 1.60217662e-19;
planck = 6.62607004e-34;
phi0 = planck/(4*pi*charge);

% resistance in kOhms

% % values from Andrei original 
% R_1 = mean([16.96]);
% R_2 = mean([14.71]);
% R_3 = mean([13.2]);
% R_4 = mean([9.83]);
% R = [R_1, R_2, R_3, R_4];
% area=[0.2, 0.25, 0.3, 0.4]*0.1; % microns^2

% % Jake FDH04 C2 values
% R = [10.30032362, 8.2, 7.30045662, 5.50060606];
% area = [0.2, 0.225, 0.3, 0.5]*0.1; % microns^2

% Jake FDH05 C1 values
R = [21.7, 16.35, 13.92, 11.6, 11.35];
area = [0.15, 0.20, 0.25, 0.30, 0.35]*0.1; % microns^2

% 
% R = [10.30032362, 8.2, 7.30045662];
% area = [0.2, 0.225, 0.3]*0.1; % microns^2

Ej = 132.6./R; %Ghz
Lj = phi0^2/planck/1e9./Ej;

for ind = 1:length(R)
    disp([num2str(area(ind)/0.1) 'x0.1 um^2, R='  num2str(R(ind)) ...
    'kOhm, Ej=' num2str(Ej(ind)) 'GHz, Lj=' num2str(Lj(ind)*1e9) 'nH']);
end
disp(newline)

areafit = linspace(0,0.1,500);
P = polyfit(area,Ej,1);
Ejfit = P(1)*areafit+P(2);

figure(1); hold on; box on
plot(area,Ej,'ro','MarkerSize',8); hold on
plot(areafit,Ejfit,'-k','LineWidth',1.5)
xlabel('Area (um^2)')
ylabel('Junction energy (GHz)')

%% Estimate Lj from junction size
L= 0.1; %um length in CAD
w = 0.22; %um width
area = w*L;
% area = w*L;
Ej = P(1)*area+P(2);
Lj =  phi0^2/planck/1e9/Ej;
% plot(area,Ej,'bo','MarkerSize',10);
disp(['LJ ' num2str(Lj/1e-9) ' nH, EJ ' num2str(Ej) ' GHz : Junction area ' num2str(L) 'um x ' num2str(w) 'um'])
%%
Ej = 8.34;
L= 0.1; %um length in CAD
area = (Ej-P(2))/P(1);
w = area/L;
Lj =  phi0^2/planck/1e9/Ej;
% plot(area,Ej,'bo','MarkerSize',10);
disp(['LJ ' num2str(Lj/1e-9) ' nH, EJ ' num2str(Ej) ' GHz : Junction area ' num2str(L) 'um x ' num2str(w) 'um'])

%% Large area inductors - 2.0um long
% Num JJ = 2* Num UCs - 1
% clc
charge = 1.60217662e-19;
planck = 6.62607004e-34;
phi0 = planck/(4*pi*charge);

% % resistance in kOhms
% R_4 = mean([43.21]);
% R_3 = mean([32.15]);
% R_2 = mean([21.79]);
% R_1 = mean([10.05]);
% R = [R_1, R_2, R_3, R_4]*2/1.5;
% numJJ = 2*[5, 10, 15, 20]-1;

% R = [15.9, 33.5, 51.8, 71.600, 88.0, 176.0]*2/1.5;
% numJJ = 2*[10, 20, 30, 40, 50, 100]-1;

% FDH05 C1 values
R = [19.6, 39.55, 59.15, 79.0, 96.0, 208.0]*2/1.5;
numJJ = 2*[10, 20, 30, 40, 50, 90]-1;


EL = 132.6./R; %Ghz
L =  phi0^2/planck/1e9./EL;

for ind = 1:length(R)
    disp(['Array ' num2str(numJJ(ind)) ' JJs(1.5x0.2 um^2), R='  num2str(R(ind)) ...
    'kOhm, EL=' num2str(EL(ind)) 'GHz, L=' num2str(L(ind)*1e9) 'nH, Lj=' num2str(L(ind)/numJJ(ind)*1e9) 'nH']);
end

disp(newline)

numJJ_fit = linspace(0,100,500);
P = polyfit(numJJ,R,1);
Rfit = P(1)*numJJ_fit+P(2);
P_L = polyfit(numJJ,L,1);
Lfit = P_L(1)*numJJ_fit+P_L(2);

figure(2); hold on; box on
plot(numJJ,R,'ro','MarkerSize',8); hold on
plot(numJJ_fit,Rfit,'-k','LineWidth',1.5)
xlabel('number JJs')
ylabel('Array resistance (kOhms)')

%% Estimate Lj from junction size
EL = .45;
R_ind = 132.6/EL;
nJJ = (R_ind-P(2))/P(1);
% plot(nJJ,R_ind,'bo','MarkerSize',10);
disp(['EL ' num2str(EL) ' GHz, number of JJs ' num2str(nJJ) ', number of UCs ' num2str((nJJ+1)/2)])
%%
nJJ = 2*(2*7-1);
Rind = P(1)*nJJ + P(2);
EL = 132.6/Rind;
L =  phi0^2/planck/1e9/EL;
Lj = L/nJJ;
% plot(nJJ,Rind,'bo','MarkerSize',10);
disp(['EL=' num2str(EL) 'GHz, L=' num2str(L*1e9) 'nH, Lj=' num2str(Lj*1e9) 'nH, numJJs=' num2str(nJJ)])

%% Inductor array 2*5-1=9X 1.75um  = 8.6kOhms 
R = mean([8.6])/9;
Ej = 132.6/R
Lj = phi0^2/planck/1e9/Ej
%% Inductor array 4X 0.8um  = 9.58kOhms 
R = mean([9.58])/4;
Ej = 132.6/R
Lj = phi0^2/planck/1e9/Ej

