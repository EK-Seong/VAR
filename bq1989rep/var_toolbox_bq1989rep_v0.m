% Replication of the bivariate VAR in Blanchard and Quah (1989, AER).
% Figure 2, pag. 662
%==========================================================================
% The VAR Toolbox 3.0 is required to run this code. To get the 
% latest version of the toolboxes visit: 
% https://github.com/ambropo/VAR-Toolbox
%==========================================================================
% Ambrogio Cesa Bianchi, November 2020
% ambrogio.cesabianchi@gmail.com


%% PRELIMINARIES
% =======================================================================
clear all; clear session; close all; clc
warning off all

addpath(genpath('C:\Users\eunkyu\git\VAR-Toolbox'))

% Load data
load FRED_QD_hw1_1.mat
X = [100*(FREDQD.gdpc1-lagmatrix(FREDQD.gdpc1,1))./FREDQD.gdpc1,FREDQD.unrate];
X = X(2:end,:);

% [xlsdata, xlstext] = xlsread('BQ1989_Data.xlsx','Sheet1');
% % Define and transform (if needed)
% X = xlsdata;
% dates = xlstext(3:end,1);
% vnames_long = xlstext(1,2:end);
% vnames = xlstext(2,2:end);
% nvar = length(vnames);
% data   = Num2NaN(xlsdata);
% % Store variables in the structure DATA
% for ii=1:length(vnames)
%     DATA.(vnames{ii}) = data(:,ii);
% end
% % Convert the first date to numeric
% year = str2double(xlstext{3,1}(1:4));
% quarter = str2double(xlstext{3,1}(6));
% % Observations
% nobs = size(data,1);


%% VAR ESTIMATION
% =======================================================================
% Define number of variables and of observations
[nobs, nvar] = size(X);
% Set deterministics for the VAR
det = 1;
% Set number of nlags
nlags = 4;
% Estimate VAR
[VAR, VARopt] = VARmodel(X,nlags,det);


%% COMPUTE IRF AND FEVD
% =======================================================================
% Set options some options for IRF calculation
VARopt.nsteps = 40;
VARopt.ident = 'long';
% VARopt.vnames = vnames_long;
VARopt.FigSize = [26,12];

% Compute IR
[IR, VAR] = VARir(VAR,VARopt);
% Compute error bands
[IRinf,IRsup,IRmed,IRbar] = VARirband(VAR,VARopt);

% Compute VD
[VD, VAR] = VARvd(VAR,VARopt);

%% REPLICATE FIGURE 1 OF BLANCHARD AND QUAH
% =======================================================================
FigSize(26,8)
% Plot supply shock
subplot(1,2,1)
plot(cumsum(IR(:,1,1)),'LineWidth',2.5,'Color',cmap(1))
hold on
plot(IR(:,2,1),'LineWidth',2.5,'Color',cmap(2))
hold on
plot(zeros(VARopt.nsteps),'--k')
title('Supply shock')
legend({'GDP Level';'Unemployment'})
% Plot demand shock
subplot(1,2,2)
plot(cumsum(IR(:,1,2)),'LineWidth',2.5,'Color',cmap(1))
hold on
plot(IR(:,2,2),'LineWidth',2.5,'Color',cmap(2))
hold on
plot(zeros(VARopt.nsteps),'-k')
title('Demand shock')
legend({'GDP Level';'Unemployment'})


FigSize(26,8)
subplot(1,2,1)
plot(VD(:,2,1),'LineWidth',2.5,'Color',cmap(2))
subplot(1,2,2)
plot(VD(:,2,2),'LineWidth',2.5,'Color',cmap(1))
legend({'GDP';'Unemployment'})

FigSize(26,8)
subplot(1,2,1)
AreaPlot(VD(:,2,1));
ylim ([0,100])
ylabel 'Contribution of Demand Shock(shaded)'
title 'GDP FEVD'
subplot(1,2,2)
AreaPlot(VD(:,2,2));
ylim ([0,100])
ylabel 'Contribution of Demand Shock(shaded)'
title 'Unemployment FEVD'


VARvdplot(VD,VARopt)

rmpath(genpath('C:\Users\eunkyu\git\VAR-Toolbox')) 