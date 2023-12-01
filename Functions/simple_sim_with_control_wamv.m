function output_MS_hist = simple_sim_with_control_wamv(t_start_MS, t_step_MS, t_end_MS, ...
        n1_cmd, n2_cmd,a1_cmd,a2_cmd, ...
        x,y,psi, ...
        u,v,r)
X = [x,y,psi,u,v,r];
U = [n1_cmd,n2_cmd,a1_cmd,a2_cmd];

nStep = (t_end_MS - t_start_MS) / t_step_MS + 1;
tSpan = linspace(t_start_MS, t_end_MS, nStep);

[t_hist_MS, state_hist_MS]= ode45(@(tODE, X) B_tugModel8(tODE,X,U),tSpan,...
    X);

data_length = length(t_hist_MS);

n1_cmd_hist = n1_cmd * ones(data_length, 1);
n2_cmd_hist = n2_cmd * ones(data_length, 1);
a1_cmd_hist = a1_cmd * ones(data_length, 1);
a2_cmd_hist = a2_cmd * ones(data_length, 1);


output_MS_hist = [t_hist_MS,...
    n1_cmd_hist, n2_cmd_hist,a1_cmd_hist, a2_cmd_hist,...
    state_hist_MS]; 

end

function X_dot = B_tugModel8(~,X, U)
    n1_cmd = U(1); %left
    n2_cmd = U(2); %right
    a1_cmd = U(3); % rad
    a2_cmd = U(4);
    
    % state variables
    x = X(1);
    y = X(2);
    psi = X(3); % rad
    u  = X(4);
    v  = X(5);
    r  = X(6); % rad
    

    
%     omega = sqrt(2300/(1000*0.004422*0.2^4));
    f1 = calThrusterForce(n1_cmd);
    f2 = calThrusterForce(n2_cmd);
    
    %注意因为输入的命令都是角度，而使用SX的时候不支持cosd，因此改为弧度。
    
    [tau_u,tau_v,tau_r] = calTorques(a1_cmd,f1,a2_cmd,f2); %这里是错的
    
    Tau = [tau_u;
        tau_v;
        tau_r];
    
    % calculate all variable dots 
    R_psi = [cos(psi), -sin(psi),0;
        sin(psi),cos(psi),0;
        0,0,1];
    
    V = [u;
        v;
        r];  
    % parameters vector
%     ref: vrx_urdf/wamv_gazebo/urdf/dynamics/wamv_gazebo_dynamics_plugin.xacro
    ship.mass = 180.0;
    ship.ixx = 120.0;
    ship.ixy = 0.0;
    ship.ixz = 0.0;
    ship.iyy = 393.0;
    ship.iyz = 0.0;
    ship.izz = 446.0;

%     vrx_urdf/wamv_gazebo/urdf/dynamics/wamv_gazebo_dynamics_plugin.xacro
    ship.xDotU = 0.0;
    ship.yDotV = 0.0;
    ship.nDotR = 0.0;

    ship.xU = -100.0;
    ship.xUU = -150.0;
    ship.yV = -100.0;
    ship.yVV = -100.0;
    ship.zW = -500.0;
    ship.kP = -300.0;
    ship.kPP = -600.0;
    ship.mQ = -900.0;
    ship.mQQ = -900.0;
    ship.nR = -800.0;
    ship.nRR = -800.0;

    ship.m11= ship.mass -ship.xDotU;
    ship.m22= ship.mass - ship.yDotV;

    ship.m33= ship.izz - ship.nDotR;

    p = ship;
    
    % dynamics
    M = [p.m11,0,0;
        0,p.m22,0;
        0,0,p.m33];
    
    c13 = -p.m22*v;
    c23 = p.m11*u;
    c31 = -c13;
    c32 = -c23;
    CV = [0,0,c13;
        0,0,c23;
        c31,c32,0];
    
    d11 = -p.xU -p.xUU*abs(u);
    d22 =-p.yV-p.yVV*abs(v);
    d33 =-p.nR-p.nRR*abs(r);
    
    DV = [d11,0,0;
        0,d22,0;
        0,0,d33];
    
    V_dot = M\(-CV*V - DV*V + Tau) ;
    
    Eta_dot = R_psi * V;
    
    % 共计13个 导数
    X_dot = [
        Eta_dot;
        V_dot];
end


%% 额外的函数，计算推力和力矩分配。
function force = calThrusterForce(rps) % 计算推力大小\
    x_u = 51.3;
    x_uu = 72.4;
    max_velocity_mps = 7.71667;
    
    max_thruster_cmd = ((x_u + x_uu * max_velocity_mps)*max_velocity_mps)/2; %2353.53
    min_thruster_cmd = -1000;
    if rps>= min_thruster_cmd && rps <= max_thruster_cmd
        force = rps;
    elseif rps > max_thruster_cmd 
        force = max_thruster_cmd;
    else
        force = min_thruster_cmd;

    end
end

function [tau_u,tau_v,tau_r] = calTorques(a1,F1,a2,F2) % 计算力矩
lx1 = -2.373776;
lx2 = -2.373776;
ly1 = 1.027135;
ly2 = -1.027135; %证明是右手坐标系。

Ta = [cos(a1), cos(a2);
    sin(a1), sin(a2);
    lx1*sin(a1)+ly1*cos(a1), lx2*sin(a2) + ly2*cos(a2)];

Tau = Ta * [F1;F2];

tau_u = Tau(1);
tau_v = Tau(2);
tau_r = Tau(3);
end