% Simulate
[t_nom, x_nom] = ode45(@(t,x) odefun(t,x,f,g,kappa_ctrl,r), tspan, x0);
t_nom = t_nom';
x_nom = x_nom';

% Compute input signal
[~,u_nom] = cellfun(@(t,x)  odefun(t,x,f,g,kappa_ctrl,r),...
    num2cell(t_nom), num2cell(x_nom,1),'UniformOutput',0);
u_nom = cell2mat(u_nom);


function [xdot,u] = odefun(t,x,f,g,kappa_ctrl,r)
  u = kappa_ctrl(x,r);
  xdot = f(x) + g(x)*u;
  disp(strcat("Simulating NOM at t = ",num2str(t)));
end