close all
clear all
clc
%%
addpath('Functions')




% Simulated source parameters
% true source
s.Q = 5; % Release rate per time unit
% source coodinates
s.x = 40; % source x coordinate
s.y = 60; % source y coordinate
s.z= 1; % hight
s.u = 4; % wind speed
s.phi = 90 * pi/180; % wind direction 
s.ci = 1; % ds 公式（1）中的危险扩散率 1
s.cii = 8; % ts 发射材料的寿命 8


% Create rectangular domain area and draw the plume 
xmin = 0;
xmax = 75;
ymin = 0;
ymax = 75;
zmin = 0;
zmax = 4;
domain = [xmin xmax ymin ymax zmin zmax]; % Size of search area

% Plot example dispersion from true source
fig1 = figure;
hold on;
hdl_plume = drawPlume(fig1, s, domain); % Draw function

% sensor model parameters 
m.thresh = 5e-4; % sensor threshold
m.Pd = 0.7; % probability of detection
m.sig = 1e-4; % minimun sensor noise
m.sig_pct = 0.5; % the standard deviation is a percentage of the concentration level

% process noise parameters (not used)
sigma.x = 0.2;
sigma.y = 0.2;
sigma.z = 0.1;
sigma.Q = 0.2;
sigma.u = 0.2;
sigma.phi = 2*pi/180;
sigma.ci = 0.1;
sigma.cii = 0.5;


% Initialisation and parameters of the mobile sensor
StartingPosition = [2 2 4]; % Starting position [x,y,z]
% change 
p.StartingState = [2,2,pi/4,0,0,0]; % x,y,/psi,u,v,r
p.StartingStateZ = 4;
p.moveT = 5;

moveDist = 2; % How far to move for one step

p.P_k = [p.StartingState, p.StartingStateZ]; % Current robot/sensor position
p.P_k_store = p.P_k;

% robot 的位置
pos.x_matrix = p.P_k(1);
pos.y_matrix = p.P_k(2);
pos.z_matrix = p.P_k(7);

D=[]; % store sensor readings


% 直接把所有点所能观测到的位置算出来，然后当运行到该点的时候能直接观测计算。
% initialise PF
N = 20000; % number of particles 

% Uniform prior for location 
theta.x = xmin + (xmax-xmin) * rand(N,1); 
theta.y = ymin + (ymax-ymin) * rand(N,1);
theta.z = zmin + (zmax-zmin) * rand(N,1);

% Gamma prior for release rate Q
a = ones(N,1)*2;
b = ones(N,1)*5;
theta.Q = gamrnd(a,b);%200*rand(N,1);%
% 这个污染源是随时间变化的gama分布。

% randn标准正太随机分布
theta.u =s.u + randn(N,1)*2;%2+6*rand(N,1);%0.75+0.5*rand(N,1);0 + randn(N,1)*0.5;
% theta.phi = s.phi + randn(N,1)*10.*pi/180;%(10 + 30*rand(N,1)).*pi/180;

% 注意phi是风的方向，先做了固定的偏移再加上随机扰动。
theta.phi = s.phi*0.9 + randn(N,1)*10.*pi/180;%(10 + 30*rand(N,1)).*pi/180;
% ci和cii还不清楚是什么。

theta.ci = s.ci+2*rand(N,1);%0.12+0.1*rand(N,1);
theta.cii = s.cii + 2*rand(N,1) - 2;%0.5+ 0.1*rand(N,1);

% 平均化固定值。
Wpnorm = ones(N,1)/N;


fig2 = figure;
hold on
preprocess(s,theta);
p.hist_ship_states = [];

for i = 1:50
    % 会找100步。
    % generate sensor data with added noise and miss-detection
    Dsim = sensorModel(s, pos, m); %污染源参数+位置+传感器参数


    D(i)=Dsim; %传感器在i步下的读数。传感器的读数其实就是决定了环境模型。

    f_dyn = @fDyn;
    h_likeli = @(s, yObv, m) hLikePlume(s, yObv, pos, m);
    g_const = @gCon; 

    [theta, Wpnorm, info] = mcmcPF(i, theta, Wpnorm, Dsim, f_dyn, sigma, h_likeli, m, g_const,[]);
    % Markov chain Monte Carlo particle filter = mcmcPF
    % 基于马尔科夫链的蒙特卡洛粒子滤波器。

    
    figure(fig1)
    hold off
    drawPlume(fig1, s, domain);
    hold on
    S(i) = 5+ceil(D(i)*5e4); % size of the marker
    scatter3(theta.x,theta.y,theta.z,3,'g','filled');
    plot3(pos.x_matrix,pos.y_matrix,pos.z_matrix,'ro','MarkerFaceColor','r','MarkerSize',5);
    plot3(p.P_k_store(:,1),p.P_k_store(:,2),p.P_k_store(:,7),'r-');
    Z_plotShipXY(p.P_k(1),p.P_k(2),p.P_k(3));
    scatter3(p.P_k_store(:,1),p.P_k_store(:,2),p.P_k_store(:,7),S,'r','MarkerFaceColor','red');
    view(0,90)
    
    drawnow



    % define the action set
    % 1 8 个控制方向 * 3
    % ynew = [[0,-moveDist,-moveDist,-moveDist,0,+moveDist,+moveDist,+moveDist] 2*[0,-moveDist,-moveDist,-moveDist,0,+moveDist,+moveDist,+moveDist] 3*[0,-moveDist,-moveDist,-moveDist,0,+moveDist,+moveDist,+moveDist]]; % [0,moveDist,-moveDist,0]
    % xnew = [[moveDist,moveDist,0,-moveDist,-moveDist,-moveDist,0,+moveDist] 2*[moveDist,moveDist,0,-moveDist,-moveDist,-moveDist,0,+moveDist] 3*[moveDist,moveDist,0,-moveDist,-moveDist,-moveDist,0,+moveDist]]; % [moveDist,0,0,-moveDist]
    % znew = [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0]; % [0,0,0,0]

    
    p.control_n1 = linspace(-300,300,11);
    p.control_n2 = linspace(-300,300,11);
    p.control_a1 = 0;
    p.control_a2 = 0;

    

    for ii = 1:size(p.control_n1,2)
        for jj = 1:size(p.control_n2,2)
            hist_ship_states = simple_sim_with_control_wamv(0,0.1,p.moveT,p.control_n1(ii),p.control_n2(jj),0,0, ...
            p.P_k(1),p.P_k(2),p.P_k(3),p.P_k(4),p.P_k(5),p.P_k(6));
            p.hist_ship_states = [p.hist_ship_states; hist_ship_states];
            p.xnew(i,ii*jj) = hist_ship_states(end,6);
            p.ynew(i,ii*jj) = hist_ship_states(end,7);
            p.znew(i,ii*jj) = 0;
            p.psinew(i,ii*jj) = hist_ship_states(end,8);
            p.unew(i,ii*jj) = hist_ship_states(end,9);
            p.vnew(i,ii*jj) = hist_ship_states(end,10);
            p.rnew(i,ii*jj) = hist_ship_states(end,11);
            
        end
    end



    % 2 4 个控制方向 
    % ynew = [[0,moveDist,-moveDist,0] 2*[0,moveDist,-moveDist,0] 3*[0,moveDist,-moveDist,0]]; % [0,moveDist,-moveDist,0]
    % xnew = [[moveDist,0,0,-moveDist] 2*[moveDist,0,0,-moveDist] 3*[moveDist,0,0,-moveDist]]; % [moveDist,0,0,-moveDist]
    % znew = [0,0,0,0, 0,0,0,0, 0,0,0,0]; % [0,0,0,0]

    % 3 只能前进，不能后退

    % 4 结合运动学模型的运动方向。

    
    % 
    Xneighbour = zeros(1,size(p.xnew,2));
    Yneighbour = zeros(1,size(p.ynew,2));
    Zneighbour = zeros(1,size(p.znew,2));

    
    
    Nz= 25; % down sample the source term particles (theta_i, i=1,...N) from N to Nz for generating the hypothetical measurements
    MM = 1; % the number of hypothetical measurements for each source term particle due to measurement noise
    
    % down sample the source term particles
    [~, indx_z]= resamplingIndex(Wpnorm,Nz);



    reward = zeros(1, size(p.xnew,2));
   
    for k = 1:size(p.xnew,2)

        Xneighbour(k) = p.xnew(k);
        Yneighbour(k) = p.ynew(k);
        Zneighbour(k) = p.znew(k);
        
        if p.xnew(k)<xmin || ...
            p.xnew(k)>xmax || ...
            p.ynew(k)<ymin || ... 
            p.ynew(k)>ymax || ... 
            p.znew(k)<zmin ||  ...
            p.znew(k)>zmax
            reward(k)=NaN;
            continue
        end


        npos.x_matrix = p.xnew(k);
        npos.y_matrix = p.ynew(k);
        npos.z_matrix = p.znew(k);

        infoGain=0;

        for jj = 1:Nz
            
            d.x = theta.x(indx_z(jj));
            d.y = theta.y(indx_z(jj));
            d.z = theta.z(indx_z(jj));
            d.Q = theta.Q(indx_z(jj));
            d.u = theta.u(indx_z(jj));
            d.phi = theta.phi(indx_z(jj));
            d.ci = theta.ci(indx_z(jj));
            d.cii = theta.cii(indx_z(jj));


            for jjj = 1:MM
                
                z = sensorModel(d, npos, m); % hypothetical measurements
                
                zUpate = hLikePlume(theta, z, npos, m);
                zWp = Wpnorm.*zUpate;
                zWpnorm = zWp./sum(zWp);
            
                WW = zWpnorm./Wpnorm;
                WW(WW<=0)=1;
                WW(isinf(WW))=1;
                WW(isnan(WW))=1;
                
                % Caculate the information gain 
                % comment/uncomment to choose one of those information gains
                % Note: here we used the sum rather than the averaged value

                % ---------------------------------------------------------
                % KLD
                infoGain = infoGain + (sum(zWpnorm.*log(WW)));
                %----------------------------------------------------------
                
                %----------------------------------------------------------
                % Entropy %基于熵的信息增益
                % infoGain = infoGain - (-sum(zWpnorm.*log2(zWpnorm+(zWpnorm==0))));
                %----------------------------------------------------------

                %----------------------------------------------------------
                % Dual control reward
                % 
                % [~, indx] = resamplingIndex(zWpnorm,round(N/5)); % downsample for quick calculation 
                % posPlus = [theta.x(indx), theta.y(indx)]';
                % posPlus_avg = mean(posPlus,2);
                % covPlus = cov(posPlus');
                
                % err_x = posPlus_avg(1)-npos.x_matrix;
                % err_y = posPlus_avg(2)-npos.y_matrix; 
                
                % infoGain = infoGain - ((err_x^2+err_y^2) + trace(covPlus));
                %----------------------------------------------------------
            end
        end

         
        reward(k) = infoGain;
        
    end


    [val,ind] = max(reward);
    
    pos.x_matrix = Xneighbour(ind);
    pos.y_matrix = Yneighbour(ind);
    pos.z_matrix = Zneighbour(ind);

    % 在这里加入船的模型，然后实时更新位置。
    % p.P_k = [pos.x_matrix pos.y_matrix pos.z_matrix];
    
    p.P_k = [p.xnew(ind),p.ynew(ind),p.psinew(ind),p.unew(ind),p.vnew(ind),p.rnew(ind),p.znew(ind)];
    p.P_k_store = [p.P_k_store; p.P_k];
   

    % stop criteria 
    [~, indx] = resamplingIndex(Wpnorm);
    Covar = cov(theta.x(indx),theta.y(indx));
    Spread = sqrt(trace(Covar));
    
    if Spread<5 % 3.5? 4?
        break
    end

end



[~, indx] = resamplingIndex(Wpnorm);

theta.x = theta.x(indx);
theta.y = theta.y(indx);
theta.z = theta.z(indx);
theta.Q = theta.Q(indx);
theta.u = theta.u(indx);
theta.phi = theta.phi(indx);
theta.ci = theta.ci(indx);
theta.cii = theta.cii(indx);

figure(fig2)
hold on
preprocess(s,theta);








