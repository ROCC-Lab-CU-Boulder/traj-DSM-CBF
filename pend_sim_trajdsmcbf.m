global isTerminalGlobal_trajdsmcbf;
setIsTerminalGlobal(false);

% Parameters
eta = 0.1;
M_slack = 1e8;
Gammabar_trajdsmcbf_sym = 400; %8
alpha1_trajdsmcbf = 100; % x_max - x >= 0
alpha2_trajdsmcbf = 100; % x_max + x >= 0
alpha3_trajdsmcbf = 100; % theta_max + pi - theta >= 0
alpha4_trajdsmcbf = 100; % theta_max - pi + theta >= 0
alpha5_trajdsmcbf = 100; % u_max - pi(x,v) >= 0
alpha6_trajdsmcbf = 100; % u_max + pi(x,v) >= 0
% Terminal DSM alphas
alpha7_trajdsmcbf = 400; % x_max - x >= 0
alpha8_trajdsmcbf = 400; % x_max + x >= 0
alpha9_trajdsmcbf = 400; % theta_max + pi - theta >= 0
alpha10_trajdsmcbf = 400; % theta_max - pi + theta >= 0
alpha11_trajdsmcbf = 400; % u_max - pi(x,v) >= 0
alpha12_trajdsmcbf = 400; % u_max + pi(x,v) >= 0
alpha13_trajdsmcbf = 400; % Gammabar


Alpha_trajdsmcbf = diag([alpha1_trajdsmcbf, alpha2_trajdsmcbf,...
    alpha3_trajdsmcbf, alpha4_trajdsmcbf, alpha5_trajdsmcbf,...
    alpha6_trajdsmcbf]);
AlphaT_trajdsmcbf = diag([alpha7_trajdsmcbf, alpha8_trajdsmcbf,...
    alpha9_trajdsmcbf, alpha10_trajdsmcbf, alpha11_trajdsmcbf,...
    alpha12_trajdsmcbf, alpha13_trajdsmcbf]);

rho_trajdsmcbf = @(v) 1*(r-v);

% Symbolic Expressions
fpi_sym = simplify(f_sym + g_sym*pi_ctrl_sym);
dfpidx_sym = simplify(jacobian(fpi_sym,xvec_sym));
dfpidv_sym = simplify(jacobian(fpi_sym,v_sym));
omega_sym = [x_max - x_sym; x_max + x_sym; ...
    theta_max + pi - theta_sym; theta_max - pi + theta_sym;...
    u_max - pi_ctrl_sym; u_max + pi_ctrl_sym];
domegadx_sym = simplify(jacobian(omega_sym,xvec_sym));
domegadv_sym = simplify(jacobian(omega_sym,v_sym));
Vfun_trajdsmcbf_sym = (xvec_sym - x_bar(v_sym))'*P_pi*(xvec_sym - x_bar(v_sym));
H_trajdsmcbf_sym = [[1,0,0,0;-1,0,0,0;0,1,0,0;0,-1,0,0]; -K_pi; K_pi];
delta_trajdsmcbf_sym = [x_max; x_max; theta_max + pi; theta_max - pi;...
    u_max - K_pi*x_bar_sym; u_max + K_pi*x_bar_sym];
Gammastar_trajdsmcbf_sym = ((delta_trajdsmcbf_sym - H_trajdsmcbf_sym*x_bar_sym).*...
    abs((delta_trajdsmcbf_sym - H_trajdsmcbf_sym*x_bar_sym)))...
    ./diag(H_trajdsmcbf_sym/P_pi*H_trajdsmcbf_sym');
DeltaT_trajdsmcbf_sym = [Gammastar_trajdsmcbf_sym; Gammabar_trajdsmcbf_sym] - Vfun_trajdsmcbf_sym;
dDeltaTdx_trajdsmcbf_sym = simplify(jacobian(DeltaT_trajdsmcbf_sym,xvec_sym));
dDeltaTdv_trajdsmcbf_sym = simplify(jacobian(DeltaT_trajdsmcbf_sym,v_sym));


% Implicit Functions
fpi = matlabFunction(fpi_sym,'vars',{xvec_sym,v_sym});
dfpidx = matlabFunction(dfpidx_sym,'vars',{xvec_sym,v_sym});
dfpidv = matlabFunction(dfpidv_sym,'vars',{xvec_sym,v_sym});
omega = matlabFunction(omega_sym,'vars',{xvec_sym,v_sym});
domegadx = matlabFunction(domegadx_sym,'vars',{xvec_sym,v_sym});
domegadv = matlabFunction(domegadv_sym,'vars',{xvec_sym,v_sym});
Vfun_trajdsmcbf = matlabFunction(Vfun_trajdsmcbf_sym,'vars',{xvec_sym,v_sym});
Gammastar_trajdsmcbf = matlabFunction(Gammastar_trajdsmcbf_sym,'vars',{v_sym});
DeltaT_trajdsmcbf = matlabFunction(DeltaT_trajdsmcbf_sym,'vars',{xvec_sym,v_sym}); 
dDeltaTdx_trajdsmcbf = matlabFunction(dDeltaTdx_trajdsmcbf_sym,'vars',{xvec_sym,v_sym}); 
dDeltaTdv_trajdsmcbf = matlabFunction(dDeltaTdv_trajdsmcbf_sym,'vars',{xvec_sym,v_sym}); 

% CBF-QP constant functions
H_qp = [2 0 0; 0 2*eta 0; 0 0 0];
f_qp = @(x,v) [-2*kappa_ctrl(x,r); -2*eta*rho_trajdsmcbf(v); M_slack];


% Simulate
Opt    = odeset('Events', @stopODE);
[t_trajdsmcbf, xv_trajdsmcbf] = ode45(@(t,xv) odefun(t,xv,fpi,...
    dfpidx,dfpidv,tspan_pred,f,g,omega,domegadx,domegadv,...
    u_max,Alpha_trajdsmcbf, H_qp, f_qp, DeltaT_trajdsmcbf,...
    dDeltaTdx_trajdsmcbf, dDeltaTdv_trajdsmcbf, AlphaT_trajdsmcbf),...
    tspan, [x0;v0], Opt);
t_trajdsmcbf = t_trajdsmcbf';
xv_trajdsmcbf = xv_trajdsmcbf';
x_trajdsmcbf = xv_trajdsmcbf(1:4,:);
v_trajdsmcbf = xv_trajdsmcbf(5,:);

% Compute input signal
[~,uws_trajdsmcbf] = cellfun(@(t,xv)  odefun(t,xv,fpi,dfpidx,dfpidv,...
    tspan_pred,f,g, omega,domegadx,domegadv,u_max,Alpha_trajdsmcbf,...
    H_qp, f_qp, DeltaT_trajdsmcbf, dDeltaTdx_trajdsmcbf, ...
    dDeltaTdv_trajdsmcbf, AlphaT_trajdsmcbf),...
    num2cell(t_trajdsmcbf), num2cell(xv_trajdsmcbf,1),...
    'UniformOutput',0);
uws_trajdsmcbf = cell2mat(uws_trajdsmcbf);
u_trajdsmcbf = uws_trajdsmcbf(1,:);
w_trajdsmcbf = uws_trajdsmcbf(2,:);
s_trajdsmcbf = uws_trajdsmcbf(3,:);

function [xvdot, uws] = odefun(t,xv,fpi,dfpidx,dfpidv,tspan_pred,f,g, ...
    omega,domegadx,domegadv,u_max,Alpha_trajdsmcbf, H_qp, f_qp, ...
    DeltaT_trajdsmcbf, dDeltaTdx_trajdsmcbf, dDeltaTdv_trajdsmcbf, ...
    AlphaT_trajdsmcbf)
  x = xv(1:4);
  v = xv(5);
  % Compute prediction trajectory
  [Sx_pred, Sv_pred, x_pred, t_pred] = ...
          dsm_compute_sensitivity_jacobians(...
          x,v,fpi,dfpidx,dfpidv,[0, tspan_pred(end)]);

  % Prepare CBF-QP (all taus)
  A_qp = zeros(6*length(t_pred),3);
  b_qp = zeros(6*length(t_pred),1);
  for i = 1:length(t_pred)
    A_qp(6*i-5:6*i,:) = -[domegadx(x_pred(:,i),v)*Sx_pred(:,:,i)*g(x), ...
        domegadx(x_pred(:,i),v)*Sv_pred(:,:,i) + domegadv(x_pred(:,i),v),...
        ones(6,1)];
    b_qp(6*i-5:6*i,:) = Alpha_trajdsmcbf*omega(x_pred(:,i), v) + ...
        domegadx(x_pred(:,i),v)*Sx_pred(:,:,i)*f(x);
  end
  % Terminal cst
  A_qp_T = -[dDeltaTdx_trajdsmcbf(x_pred(:,end),v)*Sx_pred(:,:,end)*g(x), ...
             dDeltaTdx_trajdsmcbf(x_pred(:,end),v)*Sv_pred(:,:,end) ...
           + dDeltaTdv_trajdsmcbf(x_pred(:,end),v), ...
           ones(7,1)];
  b_qp_T = AlphaT_trajdsmcbf*DeltaT_trajdsmcbf(x_pred(:,end),v) ...
         + dDeltaTdx_trajdsmcbf(x_pred(:,end),v)*Sx_pred(:,:,end)*f(x);


  % Solve CBF-QP
  [uws, ~, exitflag] = quadprog(H_qp,f_qp(x,v),[A_qp; A_qp_T],...
      [b_qp; b_qp_T],[],[],[-u_max; -inf; 0],[u_max; inf; inf]);
  % [uws, ~, exitflag] = quadprog(H_qp,f_qp(x,v), A_qp,...
  %     b_qp,[],[],[-u_max; -inf; 0],[u_max; inf; inf]);
  if exitflag == -2
    setIsTerminalGlobal(true);
  end
  disp(strcat("Simulating Traj-DSM-CBF at t = ",num2str(t)," exitflag = ",string(exitflag)));
  u = uws(1);
  w = uws(2);
  xvdot = [f(x) + g(x)*u; w];
end


function setIsTerminalGlobal(val)
  global isTerminalGlobal_trajdsmcbf
  isTerminalGlobal_trajdsmcbf = val;
end

function [value, isterminal, direction] = stopODE(t, x)
  global isTerminalGlobal_trajdsmcbf;
  value      = isTerminalGlobal_trajdsmcbf;
  isterminal = 1;   % Stop the integration if true
  direction  = 0;
end