% General parameters
r = 4; % desired reference
tspan = [0,8]; % Simulation span
tspan_pred = [0,10]; % Prediction span
x0 = [0; pi; 0; 0]; % Initial state
v0 = 0; % Initial reference state


% Constraints
x_max = 4.5;
theta_max = deg2rad(20);
u_max = 20;

% Symbolic Expressions
mc = 1; mp = 0.5; L = 0.7; grav = 9.81;
syms x_sym theta_sym dx_sym dtheta_sym v_sym real
q_sym = [x_sym; theta_sym];
dq_sym = [dx_sym; dtheta_sym];
xvec_sym = [q_sym; dq_sym];
M_sym = [mc + mp, mp*L*cos(theta_sym);
    mp*L*cos(theta_sym), mp*L^2];
C_sym = [0, -mp*L*dtheta_sym*sin(theta_sym); 0, 0];
G_sym = [0; mp*grav*L*sin(theta_sym)];
B = [1; 0];
f_sym = [dq_sym; -M_sym\(C_sym*dq_sym + G_sym)];
g_sym = [zeros(2,1); M_sym\B];
x_bar_sym = [v_sym; pi; 0; 0];
u_bar_sym = 0;
Al_sym = jacobian(f_sym,xvec_sym);
Bl_sym = g_sym;


% Implicit functions
f = matlabFunction(f_sym,'var',{xvec_sym});
g = matlabFunction(g_sym,'var',{xvec_sym});
x_bar = matlabFunction(x_bar_sym,'var',{v_sym});
Al = matlabFunction(Al_sym,'var',{xvec_sym});
Bl = matlabFunction(Bl_sym,'var',{xvec_sym});

% Controller design
[K_pi, P_pi] = lqr(double(Al(x_bar(v_sym))),double(Bl(x_bar(v_sym))),diag([1,1,1,1]),5);
% K_kappa = lqr(double(Al(x_bar(v_sym))),double(Bl(x_bar(v_sym))),diag([1,1,1,1]),0.01);
K_kappa = [-35, 150, -20, 50];
pi_ctrl_sym = -K_pi*(xvec_sym - x_bar(v_sym));
pi_ctrl = @(x,v) -K_pi*(x - x_bar(v)); % Prestabilizing controller
kappa_ctrl_sym = -K_kappa*(xvec_sym - x_bar(v_sym));
kappa_ctrl = @(x,v) -K_kappa*(x - x_bar(v));







