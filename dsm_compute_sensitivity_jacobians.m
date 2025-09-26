function [S_x, S_v, x, t] = dsm_compute_sensitivity_jacobians(x0,v,f,A,B,tspan)
%DSM_COMPUTE_SENSITIVITY_JACOBIANS computes the sensitivity jacobians Q_x,
% S_v and predicted trajectory x of a closed-loop system over tspan.
%  [S_x,S_v,x,t] = DSM_COMPUTE_SENSITIVITY_MATRIX(x0,f,A,T) Returns the time-series:
%  [S_x,t] as the state sensitivity Jacobian Q_x(t),
%  [S_v,t] as the reference sensitivity Jacobian Q_v(t), and
%  [x, t] as the trajectory prediction.
%
%  x0 is the initial state vector to compute x(t) for.
%  f is an implicit function of x and v describing the prestabilizing dynamics.
%  A is an implicit function of x and v describing the Jacobian of
%  the prestabilizing dynamics with respect to x.
%  B is an implicit function of x and v describing the Jacobian of
%  the prestabilizing dynamics with respect to v.
%  tspan is the time period [0,T] to compute Q(t) for.
%
%   Copyright (c) 2025, University of Colorado Boulder

% Get state and reference dimensions
n = length(x0);
l = length(v);
x0 = reshape(x0,n,1);
v = reshape(v, l, 1);

% Initialize flattened variable z = [x,Q_x,Q_v]
z0 = zeros(n*(n+l+1),1);
z0(1:n,1) = x0; % Initial condition is current x
z0(n+1:n*(n+1),1) = reshape(eye(n),[],1); % Q_x(0) = I_n

% Simulate flattened variable
[t, z] = ode45(@(t,z) odefun(t,z,v,f,A,B,n,l),tspan, z0);

% Reshape
z = z';
x = z(1:n,:);
S_x = zeros(n,n,length(t));
S_v = zeros(n,l,length(t));
for i = 1:length(t)
  S_x(:,:,i) = reshape(z(n+1:n*(n+1),i),n,n);
  S_v(:,:,i) = reshape(z(n*(n+1)+1:end,i),n,l);
end



end

function [zdot] = odefun(t,z,v,f,A,B,n,l)
  % disp(strcat("Computing sensitivity Jacobians Q_x(t) at t = ",num2str(t)));  
  zdot = zeros(length(z),1);
  zdot(1:n,1) = f(z(1:n,1),v); % Prestabilizing dynamics
  zdot(n+1:n*(n+1),1) = reshape(A(z(1:n,1),v)...
      *reshape(z(n+1:n*(n+1),1),n,n),[],1); % State Jacobian dynamics
  zdot(n*(n+1)+1:end,1) = reshape(reshape(A(z(1:n,1),v)...
      *reshape(z(n*(n+1)+1:end,1),n,l),[],1)...
      + B(z(1:n,1),v),[],1); % Reference Jacobian dynamics
end