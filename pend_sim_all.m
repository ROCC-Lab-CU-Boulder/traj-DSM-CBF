clear; clc; close all;

pend_data(); % Load general data

% Options
plot_cst = true;
plot_nom = true;
plot_trajdsmcbf = true;

if plot_nom
  pend_sim_nom(); % Nominal
end
if plot_trajdsmcbf
  pend_sim_trajdsmcbf(); % Trajectory DSM-CBF
end

% Plot
col_nom = [0, 0, 0]; % Black
col_trajdsmcbf = [0, 0.4470, 0.7410]; % Blue
col_cst = [1, 0, 0]; % Red

figure(1)
tiledlayout(3,1,'Padding', 'compact', 'TileSpacing', 'tight')
nexttile(1)
hold on
if plot_cst
  plot([tspan(1),tspan(end)],[x_max ,x_max],'color',col_cst)
end
plt_objects = [];
plt_labels = {};
if plot_nom
  plt_objects = [plt_objects, plot(t_nom,x_nom(1,:),'color',col_nom,'LineWidth',1.1)];
  plt_labels = [plt_labels, {"Nominal"}];
end
if plot_trajdsmcbf
  plt_objects = [plt_objects, plot(t_trajdsmcbf, x_trajdsmcbf(1,:),...
      'Color',col_trajdsmcbf,'LineWidth',2)];
  plot(t_trajdsmcbf,v_trajdsmcbf,'Color',col_trajdsmcbf,'LineStyle','--');
  plt_labels = [plt_labels, {"Traj DSM-CBF"}];
end
grid on
legend(plt_objects,plt_labels,'location','southeast')
ylabel('$x(t) ~(m)$','Interpreter','latex')
set(gca,'xticklabel',[])
ylim([-1,6])


nexttile(2)
hold on
if plot_cst
  plot([tspan(1),tspan(end)],rad2deg([pi + theta_max,pi + theta_max]),'color',col_cst)
  plot([tspan(1),tspan(end)],rad2deg([pi - theta_max,pi - theta_max]),'color',col_cst)
end
if plot_nom
  plot(t_nom,rad2deg(x_nom(2,:)),'color',col_nom,'LineWidth',1.1);
end
if plot_trajdsmcbf
  plot(t_trajdsmcbf,rad2deg(x_trajdsmcbf(2,:)),'color',col_trajdsmcbf,'LineWidth',2);
end
grid on
legend(plt_objects,plt_labels)
ylabel('$\theta(t) ~(deg)$','Interpreter','latex')
set(gca,'xticklabel',[])
ylim([140,220])


nexttile(3)
hold on
if plot_cst
  plot([tspan(1),tspan(end)],[u_max,u_max],'color',col_cst) % TODO
  plot([tspan(1),tspan(end)],-[u_max,u_max],'color',col_cst) % TODO
end
if plot_nom
  plot(t_nom,u_nom,'color',col_nom,'LineWidth',1.1);
end
if plot_trajdsmcbf
  plot(t_trajdsmcbf, u_trajdsmcbf,'color',col_trajdsmcbf,'LineWidth',2)
end
grid on
legend(plt_objects,plt_labels)
ylabel('$u(t) ~ (N)$','Interpreter','latex')
xlabel('$t~ (s)$','Interpreter','latex')
ylim([-30,30])