clear
close all
clc

addpath(genpath('C:\Users\eunkyu\git\VAR-Toolbox'))
%% Parmeter 정하기
numlag = 5;     % Number of Lag Terms
vecdim = 2;     % Dimension of Y vector
h = 39;         % horizon
with_constant = 1; % 1이면 상수항 추가 0이면 상수항 없음
with_trend = 1; % 1이면 타임트렌드 추가 0이면 없음


%% Data 불러오기


load FRED_QD_hw1_1.mat %1959~2017까지의 분기별 데이터임.
gdp = FREDQD.gdpc1;
dgdp = 100*(gdp-lagmatrix(gdp,1))./gdp;
% lgdp = log(gdp);
% dgdp = 100*lgdp-lagmatrix(lgdp,1);
unrate = FREDQD.unrate;

Y = [dgdp,unrate];
Y = Y(2:end,:);  %변화율로 변환하면서 초기값이 빠짐. 따라서 2부터 시작

T = size(Y,1);

%% SVAR estimation with Long-run Restriction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% Without trend & with constant %%%%%
if (with_trend ==0 && with_constant == 1)

    X = ones(T,numlag*vecdim+1);
    for p=1:numlag
        X(:,2*p:2*p+1)=lagmatrix(Y,p);
    end
    Y = Y(p+1:end,:); %lag term을 만들면서 p+1개의 관측치가 날아감
    X = X(p+1:end,:);
    
    numobs = size(Y,1); %관측치를 새로 정의해줌
    
    Ahat = (X'*X)\X'*Y; % 우리가 원하는 행렬의 전치행렬로 결과가 나옴
    
    resid = Y-X*Ahat;
    Sigma = resid'*resid/numobs;
    
    constant = Ahat(1,:)';
    Ahats = zeros(vecdim,vecdim,numlag);
    for p=1:numlag
        Ahats(:,:,p) = Ahat(2*p:2*p+1,:)'; %각 lag polynomial을 재배열해주기 위해 전치함
    end
    A1 = eye(vecdim); 
    for p=1:numlag
        A1 = A1-Ahats(:,:,p);
    end
    
    P = chol(A1\Sigma/A1')'; %하삼각행렬로 만들어주기 위해 전치함
    B = A1*P;
    
end

%%%%% Without constant & without trend %%%%%
if (with_trend ==0 && with_constant == 0)
    T = size(Y,1);
    X = zeros(T,numlag*vecdim);

    for p=1:numlag
    X(:,2*p-1:2*p)=lagmatrix(Y,p);
    end
    Y = Y(p+1:end,:); %lag term을 만들면서 p+1개의 관측치가 날아감
    X = X(p+1:end,:);
    numobs = size(Y,1); %관측치를 새로 정의해줌

    Ahat = (X'*X)\X'*Y; % 우리가 원하는 행렬의 전치행렬로 결과가 나옴
    
    resid = Y-X*Ahat;
    Sigma = resid'*resid/numobs;

    Ahats = zeros(vecdim,vecdim,numlag);
    for p=1:numlag
        Ahats(:,:,p) = Ahat(2*p-1:2*p,:)'; %각 lag polynomial을 재배열해주기 위해 전치함
    end

    A1 = eye(vecdim);
    for p=1:numlag
        A1 = A1-Ahats(:,:,p);
    end
    
    P = chol(A1\Sigma/A1')'; %하삼각행렬로 만들어주기 위해 전치함
    B = A1*P;

end


%%%%% With constant and trend %%%%%
if with_trend ==1
    time = 1:T;
    X = ones(T,numlag*vecdim+2);
    X(:,2) = time';
    for p=1:numlag
        X(:,2*p+1:2*p+2)=lagmatrix(Y,p);
    end
    Y = Y(p+1:end,:); %lag term을 만들면서 p+1개의 관측치가 날아감
    X = X(p+1:end,:);
    
    numobs = size(Y,1); %관측치를 새로 정의해줌
    
    Ahat = (X'*X)\X'*Y; % 우리가 원하는 행렬의 전치행렬로 결과가 나옴
    
    resid = Y-X*Ahat;
    Sigma = resid'*resid/numobs;
    
    constant = Ahat(1,:)';
    trend = Ahat(2,:)';
    Ahats = zeros(vecdim,vecdim,numlag);
    for p=1:numlag
        Ahats(:,:,p) = Ahat(2*p+1:2*p+2,:)'; %각 lag polynomial을 재배열해주기 위해 전치함
    end
    A1 = eye(vecdim);
    for p=1:numlag
        A1 = A1-Ahats(:,:,p);
    end
    
    P = chol(A1\Sigma/A1')'; %하삼각행렬로 만들어주기 위해 전치함
    B = A1*P;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% VAR(p) process를 VMA(inf) process로 표현하는 작업 %%%%%%
%%%%%%%%%%%%(Wold Decomposition) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = zeros(vecdim*numlag);
for p=1:numlag
    A(1:2,2*p-1:2*p) = Ahats(:,:,p);
end
for p=1:numlag-1
    A(2*p+1:2*p+2,2*p-1:2*p)=eye(vecdim);
end

%%% IRF to demand shock
d = [0;1];
dIRF = zeros(vecdim,h+1);
IR = zeros(vecdim*numlag,1);
IR(1:vecdim,1) = B*d;
for i=1:h+1
    dIRF(:,i) = IR(1:vecdim,1);
    IR = A*IR;
end

d_irf_gdp = dIRF(1,:);
d_irf_unrate = dIRF(2,:);

%%% IRF to supply shock
s = [1;0];
sIRF = zeros(vecdim,h+1);

IR = zeros(vecdim*numlag,1);
IR(1:vecdim,1) = B*s;
for i=1:h+1
    sIRF(:,i) = IR(1:vecdim,1);
    IR = A*IR;
end

s_irf_gdp = sIRF(1,:);
s_irf_unrate = sIRF(2,:);


%% Hansen (2022) p.521
%%% AIC 
AIC = numobs*log(det(resid'*resid/numobs))+2*vecdim*(numlag*vecdim+1);

%%
%%% IRF Graph %%%
figure(1)

subplot(1,2,1)
plot(cumsum(s_irf_gdp,2),'LineWidth',2.5,'Color','blue')
xlim ([1,h])
ylim ([-0.8,1.2])
hold on
plot(s_irf_unrate,'LineWidth',2.5,'Color','yellow')
hold on
plot(zeros(h,1),'--k')
legend({'GDP Level';'Unemployment'})
title('Supply shock')
subplot(1,2,2)
plot(cumsum(d_irf_gdp,2),'LineWidth',2.5,'Color','blue')
hold on
plot(d_irf_unrate,'LineWidth',2.5,'Color','yellow')
hold on
plot(zeros(h,1),'--k')
title('Demand shock')
legend({'GDP Level';'Unemployment'})
xlim ([1,h])
ylim ([-0.8,1.2])


% %%
% %%% Forecast Error Variance Decomposition %%%
% 
% A = zeros(vecdim*numlag);
% for p=1:numlag
%     A(1:2,2*p-1:2*p) = Ahats(:,:,p);
% end
% 
% fev_k = zeros(h+1,vecdim,vecdim); % horizon(i), shock(j), variable(k)
% fev_total = zeros(h+1,1,vecdim);
% 
% Phi = eye(vecdim*numlag);
% for i=1:h+1
%     for j = 1:vecdim
%         for k = 1:vecdim
%             fev_k(i,j,k) = (Phi(k,1:vecdim)*B(:,j))^2;
%             fev_total(i,1,k) = fev_total(i,1,k)+fev_k(i,j,k);
%         end
%     end
%     Phi=A*Phi;
% end
% 
% fevd = zeros(h+1,vecdim,vecdim);
% for j = 1:vecdim
%     for k = 1:vecdim
%         fevd(:,j,k) = cumsum(fev_k(:,j,k),1)./cumsum(fev_total(:,1,k),1);
%     end
% end
% 
% figure(2)
% FigSize(26,8)
% subplot(1,2,1)
% AreaPlot(fevd(:,2,1))
% ylabel 'Contribution of Demand Shock(shaded)'
% title 'GDP FEVD'
% hold on
% subplot(1,2,2)
% AreaPlot(fevd(:,2,2))
% ylabel 'Contribution of Demand Shock(shaded)'
% title 'Unemployment FEVD'

%%
%%% Forecast Error Variance Decomposition
%%% TA session 그대로 따라하기

d_fev_gdp = d_irf_gdp.^2;
d_fev_gdp = cumsum(d_fev_gdp,2);
s_fev_gdp = s_irf_gdp.^2;
s_fev_gdp = cumsum(s_fev_gdp,2);

d_fevd_gdp = d_fev_gdp./(d_fev_gdp+s_fev_gdp);

d_fev_unrate = d_irf_unrate.^2;
d_fev_unrate = cumsum(d_fev_unrate,2);
s_fev_unrate = s_irf_unrate.^2;
s_fev_unrate = cumsum(s_fev_unrate,2);

d_fevd_unrate = d_fev_unrate./(d_fev_unrate+s_fev_unrate);

figure(2)

subplot(1,2,1)
plot(d_fevd_gdp,'LineWidth',2.5,'Color','blue')
% AreaPlot(d_fevd_gdp)
ylabel 'Contribution of Demand Shock(shaded)'
title 'GDP FEVD'
ylim ([0,1])
xlim ([1,h])
hold on
subplot(1,2,2)
plot(d_fevd_unrate,'LineWidth',2.5,'Color','yellow')
% AreaPlot(d_fevd_unrate)
ylabel 'Contribution of Demand Shock(shaded)'
title 'Unemployment FEVD'
ylim ([0,1])
xlim ([1,h])


%%
%%% Bootstrap CIs
%%% Recursive Bootstrap

rep=1000;

if (with_trend ==0 && with_constant == 1)
    boot_save_Ahat = zeros(vecdim*numlag+1,vecdim,rep);
end

if (with_trend ==0 && with_constant == 0)
    boot_save_Ahat = zeros(vecdim*numlag,vecdim,rep);
end

if with_trend==1
    boot_save_Ahat = zeros(vecdim*numlag+2,vecdim,rep);
end

boot_save_dIRF = zeros(vecdim,h+1,rep);
boot_save_sIRF = zeros(vecdim,h+1,rep);

for r=1:rep
    % random draw from empirical distribution of centered residuals
    bootsample=randi(numobs,numobs,1); 
    resid_centered = resid-mean(resid,1);
    Ub = resid_centered(bootsample,:);
    % Constructing bootstrap sample
    Yb=X*Ahat+Ub;
    Xb=X;
    if (with_trend ==0 && with_constant == 1)
        for p=1:numlag
            Xb(:,2*p:2*p+1)=lagmatrix(Yb,p);
        end
        Xb(:,1) = ones(numobs,1);   % constant term
        for p=1:numlag
            Xb(1:(2*p-1)*p,2*p:2*p+1)=X(1:(2*p-1)*p,2*p:2*p+1);
        end
    end

    if (with_trend ==0 && with_constant == 0)
        for p=1:numlag
            Xb(:,2*p-1:2*p)=lagmatrix(Yb,p);
        end
        for p=1:numlag
            Xb(1:(2*p-1)*p,2*p-1:2*p)=X(1:(2*p-1)*p,2*p-1:2*p);
        end
    end

    if with_trend==1
        for p=1:numlag
            Xb(:,2*p+1:2*p+2)=lagmatrix(Yb,p);
        end
        Xb(:,1) = ones(numobs,1);   % constant term
        Xb(:,2) = (1:1:numobs)';    % time trend
        for p=1:numlag
            Xb(1:(2*p-1)*p,2*p+1:2*p+2)=X(1:(2*p-1)*p,2*p+1:2*p+2);
        end
    end
    
    % Saving estimates for each repetition
    Ahat_b = (Xb'*Xb)\Xb'*Yb;
    boot_save_Ahat(:,:,r)=Ahat_b;

    resid_b = Yb-Xb*Ahat_b;
    Sigma_b = resid_b'*resid_b/numobs;

    Ahats_b = zeros(vecdim,vecdim,numlag);
    if (with_trend ==0 && with_constant == 1)
        constant_b = Ahat_b(1,:)';
        for p=1:numlag
            Ahats_b(:,:,p) = Ahat_b(2*p:2*p+1,:)'; %각 lag polynomial을 재배열해주기 위해 전치함
        end 
    end

    if (with_trend ==0 && with_constant == 0)
        for p=1:numlag
            Ahats_b(:,:,p) = Ahat_b(2*p-1:2*p,:)'; %각 lag polynomial을 재배열해주기 위해 전치함
        end 
    end

    if with_trend==1
        constant_b = Ahat_b(1,:)';
        trend_b = Ahat_b(2,:)';
        for p=1:numlag
            Ahats_b(:,:,p) = Ahat_b(2*p+1:2*p+2,:)'; %각 lag polynomial을 재배열해주기 위해 전치함
        end 
    end

    A1_b = eye(vecdim);
    for p=1:numlag
        A1_b = A1_b-Ahats_b(:,:,p);
    end
    
    P_b = chol(A1_b\Sigma_b/A1_b')'; %하삼각행렬로 만들어주기 위해 전치함
    B_b = A1_b*P_b;

    %%%%%% 여기서부터는 동일
    
    A_b = zeros(vecdim*numlag);
    for p=1:numlag
        A_b(1:2,2*p-1:2*p) = Ahats_b(:,:,p);
    end
    for p=1:numlag-1
        A_b(2*p+1:2*p+2,2*p-1:2*p)=eye(vecdim);
    end
    
    %%% IRF to demand shock
    d = [0;1];
    dIRF_b = zeros(vecdim,h+1);
    IR = zeros(vecdim*numlag,1);
    IR(1:vecdim,1) = B_b*d;
    for i=1:h+1
        dIRF_b(:,i) = IR(1:vecdim,1);
        IR = A_b*IR;
    end
    boot_save_dIRF(:,:,r)=dIRF_b;
   
    d_irf_gdp_b = dIRF_b(1,:);
    d_irf_unrate_b = dIRF_b(2,:);

    %%% IRF to supply shock
    s = [1;0];
    sIRF_b = zeros(vecdim,h+1);

    IR = zeros(vecdim*numlag,1);
    IR(1:vecdim,1) = B_b*s;
    for i=1:h+1
        sIRF_b(:,i) = IR(1:vecdim,1);
        IR = A_b*IR;
    end
    boot_save_sIRF(:,:,r)=sIRF_b;

    s_irf_gdp_b = sIRF_b(1,:);
    s_irf_unrate_b = sIRF_b(2,:);

end

up95_sIRF = zeros(vecdim,h+1);
lo95_sIRF = zeros(vecdim,h+1);
up95_dIRF = zeros(vecdim,h+1);
lo95_dIRF = zeros(vecdim,h+1);

boot_save_sIRF_cumsum = cumsum(boot_save_sIRF,2);
boot_save_dIRF_cumsum = cumsum(boot_save_dIRF,2);
for i=1:h+1
    up95_sIRF(2,i)=quantile(boot_save_sIRF(2,i,:),0.975,3);
    lo95_sIRF(2,i)=quantile(boot_save_sIRF(2,i,:),0.025,3);
    up95_dIRF(2,i)=quantile(boot_save_dIRF(2,i,:),0.975,3);
    lo95_dIRF(2,i)=quantile(boot_save_dIRF(2,i,:),0.025,3);

    up95_sIRF(1,i)=quantile(boot_save_sIRF_cumsum(1,i,:),0.975,3);
    lo95_sIRF(1,i)=quantile(boot_save_sIRF_cumsum(1,i,:),0.025,3);
    up95_dIRF(1,i)=quantile(boot_save_dIRF_cumsum(1,i,:),0.975,3);
    lo95_dIRF(1,i)=quantile(boot_save_dIRF_cumsum(1,i,:),0.025,3);
end

figure(3)
% title 'Bootstrap Confidence Bands'
subplot(2,2,4)

plot(lo95_dIRF(2,:))
hold on
plot(up95_dIRF(2,:))
hold on
plot(d_irf_unrate,'LineWidth',2.5,'Color','yellow')
hold on
plot(zeros(h+1,1),'--k')
ylabel 'Unemployment'
title 'SIRF-Demand Shock'

subplot(2,2,2)

plot(lo95_dIRF(1,:))
hold on
plot(up95_dIRF(1,:))
hold on
plot(cumsum(d_irf_gdp,2),'LineWidth',2.5,'Color','yellow')
hold on
plot(zeros(h+1,1),'--k')
ylabel 'GDP level'
title 'SIRF-Demand Shock'

subplot(2,2,3)

plot(lo95_sIRF(2,:))
hold on
plot(up95_sIRF(2,:))
hold on
plot(s_irf_unrate,'LineWidth',2.5,'Color','yellow')
hold on
plot(zeros(h+1,1),'--k')
ylabel 'Unemployment'
title 'SIRF-Supply Shock'

subplot(2,2,1)

plot(lo95_sIRF(1,:))
hold on
plot(up95_sIRF(1,:))
hold on
plot(cumsum(s_irf_gdp,2),'LineWidth',2.5,'Color','yellow')
hold on
plot(zeros(h+1,1),'--k')
ylabel 'GDP level'
title 'SIRF-Supply Shock'


figure(4)

subplot(2,2,1)

plot(lo95_sIRF(1,:))
hold on
plot(up95_sIRF(1,:))
for r=1:rep
    plot(cumsum(boot_save_sIRF(1,:,r),2),"yellow")
    hold on
end
plot(median(cumsum(boot_save_sIRF(1,:,:),2),3),'LineWidth',2.5,'Color','r')
hold on
plot(zeros(h+1,1),'--k')
ylabel 'GDP level'
title 'SIRF-Supply Shock'

subplot(2,2,3)

plot(lo95_sIRF(2,:))
hold on
plot(up95_sIRF(2,:))
for r=1:rep
    plot(boot_save_sIRF(2,:,r),"yellow")
    hold on
end
plot(median(boot_save_sIRF(2,:,:),3),'LineWidth',2.5,'Color','r')
hold on
plot(zeros(h+1,1),'--k')
ylabel 'Unemployment'
title 'SIRF-Supply Shock'

subplot(2,2,2)

plot(lo95_dIRF(1,:))
hold on
plot(up95_dIRF(1,:))
for r=1:rep
    plot(cumsum(boot_save_dIRF(1,:,r),2),"yellow")
    hold on
end
plot(median(cumsum(boot_save_dIRF(1,:,:),2),3),'LineWidth',2.5,'Color','r')
hold on
plot(zeros(h+1,1),'--k')
ylabel 'GDP level'
title 'SIRF-Demand Shock'

subplot(2,2,4)

plot(lo95_dIRF(2,:))
hold on
plot(up95_dIRF(2,:))
for r=1:rep
    plot(boot_save_dIRF(2,:,r),"yellow")
    hold on
end
plot(median(boot_save_dIRF(2,:,:),3),'LineWidth',2.5,'Color','r')
hold on
plot(zeros(h+1,1),'--k')
ylabel 'Unemployment'
title 'SIRF-Demand Shock'

subplot(2,2,2)

plot(lo95_dIRF(1,:))
hold on
plot(up95_dIRF(1,:))

subplot(2,2,4)

plot(lo95_dIRF(2,:))
hold on
plot(up95_dIRF(2,:))

subplot(2,2,3)

plot(lo95_sIRF(2,:))
hold on
plot(up95_sIRF(2,:))

subplot(2,2,1)

plot(lo95_sIRF(1,:))
hold on
plot(up95_sIRF(1,:))




%%%%% Without trend & with constant %%%%%
if (with_trend ==0 && with_constant == 1)
end

%%%%% Without constant & without trend %%%%%
if (with_trend ==0 && with_constant == 0)
end



rmpath(genpath('C:\Users\eunkyu\git\VAR-Toolbox')) 
















% recursive bootstrap하다가 포기함
% rep = 1000;
% init = [dgdp(2:numlag+1,1),unrate(2:numlag+1,1)];
% init_cond = zeros(1,vecdim*numlag);
% for p=1:numlag
%     init_cond(:,2*p-1:2*p)=init(p,:);
% end
% Ahat_boot = zeros(2+vecdim*numlag,vecdim,rep);
% 
% 
% if with_trend ==1
%     init_cond = [1,0,init_cond];
%     for r=1:rep
%         bootsample = randi(numobs,numobs,1);
%         error_boot = resid(bootsample,:);
%         Yb = zeros(numobs,vecdim);
%         Yb(1,:) = init_cond*Ahat;
%         for i=2:numobs
%             Yb(i,:)=
%         end
%     end
% end


























