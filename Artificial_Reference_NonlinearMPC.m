%{
Simulation for the paper: "MPC for tracking of constrained nonlinear
                           systems" - A. Ferrramosca
Reference: "Nonlinear MPC for tracking piece-wise constant reference 
            signals" - D. Limon
%}

close all;
clear all;
clc;

%% System parameters
q = 100;            % [l/min]
Tf = 350;           % [K]
V = 100;            % [l]
rho =  1000;        % [g/l]
Cp = 0.239;         % [J/g K]
DeltaH = -5*10^4;   % [J/mol]
E_R = 8750;         % [K]
k0 = 7.2*10^10;     % [min^-1]
UA = 5*10^4;        % [J/min K]
CAf = 1;            % [mol/l]

%% System constraints
x1_min = 0; x1_max = 1;     % [mol/l];  x1 = CA
x2_min = 280; x2_max = 370; % [K];      x2 = T
x_min = [x1_min; x2_min];
x_max = [x1_max; x2_max];

u_min = 280; u_max = 370;   % [K];      u = Tc

nx = 2; nu = 1; ny = 1;

%% Simulation parameters
Ts = 0.03; % [min] % sampling time

% initial conditions
x0 = [0.2057; 370];
u0 = 300.1261;

%% MPC parameters
N = 5; % prediction/control horizon

Q = diag([1, 1/100]);
R =  1/100;

% offset-cost weight
T = 100;

%%%%------------------------------%%%%------------------------------%%%%

%% Steady state and control computation
%{
This part computes the steady state and control that corresponds to a
desired output
%}
yt = 336.9; % = x2 = T [K]      % desired output

options1 = optimoptions('fsolve', 'Display', 'iter');
z0 = [0; yt; 300];
z_sol = fsolve(@(z) steadyStateFcn(z, yt, Ts), z0, options1);

xt = z_sol(1:nx); % steady state
ut = z_sol(nx+1:end); % steady control

%%%%------------------------------%%%%------------------------------%%%%

%% Domain of Attraction computation - 1 (Tracking)
%{
This part computes the Domain of Attraction for the set of admissible
target reference.

Reference: Phase I algorithm - Section 11.4 - Convex optimization - 
           S. Boyd, L. Vandenberghe

Accuracy depends on the finest of gridding (number of testing points)
and optimization algorithm: sqp / interior point / ...

Computation burden depends on the prediction horizon
%}

% % Gridding the state space
Nsamples1 = 75;
Nsamples2 = 50;
x1grid = linspace(x1_min, x1_max, Nsamples1);
x2grid = linspace(x2_min, x2_max, Nsamples2);

% % Admissible initial states
x0Ad = [];

% % fmincon options
options = optimoptions("fmincon","Algorithm","sqp","MaxIterations",100,"Display","none");

for i = 1:Nsamples1
    for j = 1:Nsamples2
        % Display iteration
        fprintf("Iteration: %d \n", (i-1)*Nsamples2+j);
        % testing point
        xtest = [x1grid(i); x2grid(j)];
        
        % Initialize decision variable
        x_seqt0 = [reshape(repmat(xtest, 1, N+1), [], 1); 
                  reshape(repmat(u0, 1, N), [], 1);
                  [xtest; u0];
                  zeros(2*nx*(N+1), 1)];
        % Solve for optimal control
        [~, fvalt, ~, ~] = fmincon(...
            @(xi_seq) getTrackingCostDOA(xi_seq, N),... % fun
            x_seqt0,...                                                      % x0
            [], [], [], [], [], [], ...                                     % A, b, Aeq, beq, lb, ub
            @(xi_seq) getTrackingConstraintsDOA(xi_seq, xtest, Ts, N, x_max, x_min, u_max, u_min),... % nonlcon
            options);                                                       % options
        % Accept the test point if the optimal cost value is 0
        if fvalt == 0
            x0Ad = [x0Ad, xtest];
        end
    end
end

% % Create the Polyhedron from the obtained set of points
set_x0Ad = Polyhedron('V', x0Ad');
set_x0Ad.minVRep(); set_x0Ad.minHRep();
figure; hold on; box on; grid on;
plot(x0Ad(1,:), x0Ad(2,:), '*');
plot(set_x0Ad, 'color', 'red', 'alpha', 0.1);
xlabel("$C_\mathrm{A}$","Interpreter","latex");
ylabel("$T$","Interpreter","latex");
hold off;

%%%%------------------------------%%%%------------------------------%%%%
%% Domain of Attraction computation - 2 (Regulation to a single target)
%{
This part computes the Domain of Attraction for a SINGLE target 
reference only.

Reference: Phase I algorithm - Section 11.4 - Convex optimization - 
           S. Boyd, L. Vandenberghe

The objective is to compare with the previous part to
see the enlargement of Domain of Attraction.
%}

% % Admissible initial states
x0Ad2 = [];

for i = 1:Nsamples1
    for j = 1:Nsamples2
        % Display iteration
        fprintf("Iteration: %d \n", (i-1)*Nsamples2+j);
        % testing point
        xtest = [x1grid(i); x2grid(j)];
        
        % Initialize decision variable
        x_seqt0 = [reshape(repmat(xtest, 1, N+1), [], 1); 
                  reshape(repmat(ut, 1, N), [], 1);
                  zeros(2*nx*(N+1), 1)];
        % Solve for optimal control
        [~, fvalt, ~, ~] = fmincon(...
            @(xi_seq) getRegulationCostDOA(xi_seq, N),... % fun
            x_seqt0,...                                                      % x0
            [], [], [], [], [], [], ...                                     % A, b, Aeq, beq, lb, ub
            @(xi_seq) getRegulationConstraintsDOA(xi_seq, xtest, xt, Ts, N, x_max, x_min, u_max, u_min),... % nonlcon
            options);                                                       % options
        % Accept the test point if the optimal cost value is 0
        if fvalt == 0
            x0Ad2 = [x0Ad2, xtest];
        end
    end
end

% % Create the Polyhedron from the obtained set of points
set_x0Ad2 = Polyhedron('V', x0Ad2');
set_x0Ad2.minVRep(); set_x0Ad2.minHRep();

% Compare the two computed Domain of Attraction

figure; hold on; box on; grid on;
plot(xt(1), xt(2), 'o');
plot(set_x0Ad, 'color', [0.8, 0.8, 0.8], 'alpha', 0.05, 'linestyle','-', 'edgecolor', 'blue');
plot(set_x0Ad2, 'color', [0.8, 0.8, 0.8], 'alpha', 0.05, 'linestyle','-', 'edgecolor', 'red');
xlabel("$C_\mathrm{A}$","Interpreter","latex");
ylabel("$T$","Interpreter","latex");
hold off;

%%%%------------------------------%%%%------------------------------%%%%

%% Closed-loop simulation

%{
This part demonstrates the case using equality terminal constraint with
prediction horizon == control horizon.
%}

% Simulation parameters
Tsim = 45; % [s]
Nsim = Tsim/Ts; % Maximum simulation iterations

% Initialization of decision variable = [xseq; useq; xs; us]
x_seq = zeros(nx, N+1); % state sequence
x_seq(:,1) = x0;
u_seq = zeros(nu, N); % control sequence
zs = [xt; ut]; % zeros(nx+nu, 1); % artificial parameter

% Initialization of simulation state, control, target reference, and artificial reference
xSim_seq = nan(nx, Nsim); xSim_seq(:,1) = x0;
uSim_seq = nan(nu, Nsim);
yref_seq = nan(ny, Nsim);
xs_seq = nan(nx, Nsim);

% % fmincon options
options = optimoptions("fmincon","Algorithm","interior-point","MaxIterations",500,"Display","iter-detailed");

% Main simulation loop
for k = 1 : Nsim
    % Update current states:
    x0 = xSim_seq(:,k);

    % Update reference
    if k <= 500
        yt = 310;
    elseif k <= 1000
        yt = 370;
    else
        yt = 290;
    end
    yref_seq(:,k) = yt;

    % Initialize decision variable for warm-start
    x_seq0 = [reshape(x_seq, nx*(N+1), 1); 
              reshape(u_seq, nu*N, 1);
              zs];

    % Solve for optimal control
    [xi_seq_opt, fval, exitflag, output] = fmincon(...
        @(xi_seq) getCost(xi_seq, yt, N, Q, R, T),... % fun
        x_seq0,...                                                      % x0
        [], [], [], [], [], [], ...                                     % A, b, Aeq, beq, lb, ub
        @(xi_seq) getConstraints(xi_seq, x0, Ts, N, x_max, x_min, u_max, u_min),... % nonlcon
        options);                                                       % options

    % Regenerate state sequence from the optimal decision variable
    x_seqtemp = reshape(xi_seq_opt(1:nx*(N+1)), nx, N+1);
    
    % Regenerate control sequence from the optimal decision variable
    u_seq = reshape(xi_seq_opt(nx*(N+1)+1:nx*(N+1)+nu*N), nu, N);
    
    % Regenerate artificial reference from the optimal decision variable
    zs = xi_seq_opt(nx*(N+1)+nu*N+1:end);
    xs = zs(1:nx);
    us = zs(nx+1:end);

    % Store optimal control
    uSim_seq(:,k) = u_seq(:,1);
    % Store artificial reference
    xs_seq(:,k) = xs;

    % Shift solution 1 step for next iteration (warm-start)
    x_seq = [x_seqtemp(:,2:N+1), xs];
    u_seq = [u_seq(:,2:N), us];

    % Propagate system
    xSim_seq(:,k+1) = CSTR_ddyns(Ts, xSim_seq(:,k), uSim_seq(:,k));

end

%%%%------------------------------%%%%------------------------------%%%%

%% Plots
t_seq = 0:1:Nsim;
figure;
tiledlayout(2,1,"TileSpacing","tight","Padding","tight");

nexttile;
grid on; hold on; box on;
plot(t_seq, xSim_seq(1,:),'-k','LineWidth',1);
plot(t_seq, [nan,xs_seq(1,:)],'--b','LineWidth',1);
legend({"Actual","Artificial"},"Interpreter","latex");
ylim([0.1, 1]);
xlabel("Samples","Interpreter","latex");
ylabel("$C_\mathrm{A}$","Interpreter","latex");
hold off;

nexttile;
grid on; hold on; box on;
plot(t_seq, xSim_seq(2,:),'-k','LineWidth',1);
plot(t_seq, [nan,yref_seq],'--k','LineWidth',1);
plot(t_seq, [nan,xs_seq(2,:)],'-.r','LineWidth',1);
legend({"Actual","Target","Artificial"},"Interpreter","latex");
ylim([280 380]);
xlabel("Samples","Interpreter","latex");
ylabel("$T$","Interpreter","latex");
hold off;

%%%%------------------------------%%%%------------------------------%%%%

%% Function helper: Cost function for MPC design

function V_N = getCost(xi_seq, yt, N, Q, R, T)

% Dimensions
nx = 2; nu = 1;

% Regenerate state sequence from the decision variable
x_seq = reshape(xi_seq(1:nx*(N+1)), nx, N+1);

% Regenerate control sequence from the decision variable
u_seq = reshape(xi_seq(nx*(N+1)+1:nx*(N+1)+nu*N), nu, N);

% Regenerate artificial reference from the decision variable
zs = xi_seq(nx*(N+1)+nu*N+1:end);
xs = zs(1:nx);
us = zs(nx+1:end);

% Initialize cost value
V_N = 0;
% Compute cost over prediction horizon
for k = 1 : N
    % Term 1: state tracking
    term1 = (x_seq(:,k)-xs)'*Q*(x_seq(:,k)-xs);
    % Term 2: control tracking
    term2 = (u_seq(:,k)-us)'*R*(u_seq(:,k)-us);
    V_N = V_N + term1 + term2;
end
% Term 3: artificial steady state deviation from target reference
term3 = (xs(2)-yt)'*T*(xs(2)-yt);
% Total cost
V_N = V_N + term3;

end

%%%%------------------------------%%%%------------------------------%%%%

%% Function helper: Constraint function for MPC design

function [c_ineq, c_eq] = getConstraints(xi_seq, x0, Ts, N, x_max, x_min, u_max, u_min)
% Dimensions
nx = 2; nu = 1;

% Regenerate state sequence from the decision variable
x_seq = reshape(xi_seq(1:nx*(N+1)), nx, N+1);
% Regenerate control sequence from the decision variable
u_seq = reshape(xi_seq(nx*(N+1)+1:nx*(N+1)+nu*N), nu, N);
% Regenerate artificial reference from the decision variable
zs = xi_seq(nx*(N+1)+nu*N+1:end);
xs = zs(1:nx);
us = zs(nx+1:end);

% State continuity constraints
c_eq_sc = nan(nx*(N+1),1);
c_eq_sc(1:nx,1) = x_seq(:,1)-x0; 
for k = 1 : N
    xkp1 = CSTR_ddyns(Ts, x_seq(:,k), u_seq(:,k));
    c_eq_sc(nx*k+1:nx*k+nx,1) = x_seq(:,k+1)-xkp1;
end

% State hard constraints from 0 to N-1
c_ineq_x = nan(2*nx*N,1);
for k = 1 : N
    c_ineq_x(2*nx*(k-1)+1 : 2*nx*(k-1)+2*nx) = [x_min-x_seq(:,k);
                                                x_seq(:,k)-x_max];
end

% Control hard constraints from 0 to N-1
c_ineq_u = nan(2*nu*N, 1);
for k = 1 : N
    c_ineq_u(2*nu*(k-1)+1 : 2*nu*(k-1)+2*nu) = [u_min-u_seq(:,k);
                                                u_seq(:,k)-u_max];
end

% Terminal constraint at N:
c_eq_xN = x_seq(:,N+1) - xs;

% Steady-state constraints
c_eq_ss = CSTR_ddyns(Ts, xs, us) - xs;

% % Concatenate constraints
% Inequality constraints
c_ineq = [c_ineq_x; % hard constraint on state(from 0 to N-1)
          c_ineq_u]; % hard constraint on control(from 0 to N-1)
          
% Equality constraint
c_eq = [c_eq_sc; % Continuity constraints (Multiple shooting method)
        c_eq_xN;
        c_eq_ss];

end

%%%%------------------------------%%%%------------------------------%%%%

%% Function helper: Steady state and control function

function F = steadyStateFcn(z, yt, Ts)
nx = 2;

x = z(1:nx);
u = z(nx+1:end);

fx = CSTR_ddyns(Ts, x, u);
hx = x(2);

% Steady-state equations
F1 = x - fx;
F2 = hx - yt;
F = [F1; F2];

end

%%%%------------------------------%%%%------------------------------%%%%

%% Function helper: Discrete-time dynamics using 5th-order Runge Kutta

function xkp1 = CSTR_ddyns(h, x, u)

k1 = h*CSTR_cdyns(x, u);
k2 = h*CSTR_cdyns(x + 1/4*k1, u);
k3 = h*CSTR_cdyns(x + 3/32*k1 + 9/32*k2, u);
k4 = h*CSTR_cdyns(x + 1932/2197*k1 - 7200/2197*k2 + 7296/2197*k3, u);
k5 = h*CSTR_cdyns(x + 439/216*k1 - 8*k2 + 3680/513*k3 - 845/4104*k4, u);
k6 = h*CSTR_cdyns(x - 8/27*k1 + 2*k2 - 3544/2565*k3 + 1859/4104*k4 - 11/40*k5, u);

xkp1 = x + 16/135*k1 + 6656/12825*k3 + 28561/56430*k4 - 9/50*k5 + 2/55*k6;

end

%%%%------------------------------%%%%------------------------------%%%%

%% Function helper: Continuous-time dynamics
function dxdt = CSTR_cdyns(x, u)
% % System parameters
q = 100;            % [l/min]
Tf = 350;           % [K]
V = 100;            % [l]
rho =  1000;        % [g/l]
Cp = 0.239;         % [J/g K]
DeltaH = -5*10^4;   % [J/mol]
E_R = 8750;         % [K]
k0 = 7.2*10^10;     % [min^-1]
UA = 5*10^4;        % [J/min K]
CAf = 1;            % [mol/l]

% % System states
CA = x(1);
T = x(2);
% % System control input
Tc = u;

% % continuous-time dynamics
dxdt = zeros(2,1);
dxdt(1) = (q/V)*(CAf-CA)-k0*exp((-E_R/T))*CA;
dxdt(2) = (q/V)*(Tf-T)-DeltaH/(rho*Cp)*k0*exp(-E_R/T)*CA + UA/(V*rho*Cp)*(Tc-T);

end

%%%%------------------------------%%%%------------------------------%%%%

%% Function helper: Constraint function for DOA - 1

function [c_ineq, c_eq] = getTrackingConstraintsDOA(xi_seq, x0, Ts, N, x_max, x_min, u_max, u_min)
% Dimensions
nx = 2; nu = 1;

% Regenerate state sequence from the decision variable
x_seq = reshape(xi_seq(1:nx*(N+1)), nx, N+1);
% Regenerate control sequence from the decision variable
u_seq = reshape(xi_seq(nx*(N+1)+1:nx*(N+1)+nu*N), nu, N);
% Regenerate artificial reference from the decision variable
zs = xi_seq(nx*(N+1)+nu*N+1:nx*(N+1)+nu*N+nx+nu);
xs = zs(1:nx);
us = zs(nx+1:end);
% Regenerate slack variables from the decision variables
s_seq = reshape(xi_seq(nx*(N+1)+nu*N+nx+nu+1:end), 2*nx, N+1);

% State continuity constraints
c_eq_sc = nan(nx*(N+1),1);
c_eq_sc(1:nx,1) = x_seq(:,1)-x0; 
for k = 1 : N
    xkp1 = CSTR_ddyns(Ts, x_seq(:,k), u_seq(:,k));
    c_eq_sc(nx*k+1:nx*k+nx,1) = x_seq(:,k+1)-xkp1;
end

% State hard constraints from 0 to N
c_ineq_x = nan(2*nx*(N+1),1);
for k = 1 : N+1
    c_ineq_x(2*nx*(k-1)+1 : 2*nx*(k-1)+2*nx) = [x_min-x_seq(:,k);
                                                x_seq(:,k)-x_max ] - s_seq(:,k);
end

% Control hard constraints from 0 to N-1
c_ineq_u = nan(2*nu*N, 1);
for k = 1 : N
    c_ineq_u(2*nu*(k-1)+1 : 2*nu*(k-1)+2*nu) = [u_min-u_seq(:,k);
                                                u_seq(:,k)-u_max];
end

% Terminal constraint at N:
c_eq_xN = x_seq(:,N+1) - xs;

% Steady-state constraints
c_eq_ss = CSTR_ddyns(Ts, xs, us) - xs;

% Constraints on slack variables
c_ineq_s = nan(2*nx*(N+1),1);
for k = 1 : N+1
    c_ineq_s(2*nx*(k-1)+1 : 2*nx*(k-1)+2*nx) = -s_seq(:,k);
end

% % Concatenate constraints
% Inequality constraints
c_ineq = [c_ineq_x; % boundary constraint on state (from 0 to N)
          c_ineq_u; % boundary constraint on control (from 0 to N-1)
          c_ineq_s]; % non-negative constraint on slacks
          
% Equality constraint
c_eq = [c_eq_sc; % Continuity constraints (Multiple shooting method)
        c_eq_xN; % Terminal constraints: xN == xs
        c_eq_ss]; % Steady-state constraint

end

%%%%------------------------------%%%%------------------------------%%%%

%% Function helper: Cost function for DOA - 1

function JN = getTrackingCostDOA(xi_seq, N)
% Dimensions
nx = 2; nu = 1;

% Regenerate slack variables from the decision variables
s_seq = xi_seq(nx*(N+1)+nu*N+nx+nu+1:end);

% Compute cost
JN = s_seq'*s_seq;

end

%%%%------------------------------%%%%------------------------------%%%%

%% Function helper: Constraint function for DOA - 2

function [c_ineq, c_eq] = getRegulationConstraintsDOA(xi_seq, x0, xt, Ts, N, x_max, x_min, u_max, u_min)
% Dimensions
nx = 2; nu = 1;

% Regenerate state sequence from the decision variable
x_seq = reshape(xi_seq(1:nx*(N+1)), nx, N+1);
% Regenerate control sequence from the decision variable
u_seq = reshape(xi_seq(nx*(N+1)+1:nx*(N+1)+nu*N), nu, N);

% Regenerate slack variables from the decision variables
s_seq = reshape(xi_seq(nx*(N+1)+nu*N+1:end), 2*nx, N+1);

% State continuity constraints
c_eq_sc = nan(nx*(N+1),1);
c_eq_sc(1:nx,1) = x_seq(:,1)-x0; 
for k = 1 : N
    xkp1 = CSTR_ddyns(Ts, x_seq(:,k), u_seq(:,k));
    c_eq_sc(nx*k+1:nx*k+nx,1) = x_seq(:,k+1)-xkp1;
end

% State hard constraints from 0 to N
c_ineq_x = nan(2*nx*(N+1),1);
for k = 1 : N+1
    c_ineq_x(2*nx*(k-1)+1 : 2*nx*(k-1)+2*nx) = [x_min-x_seq(:,k);
                                                x_seq(:,k)-x_max ] - s_seq(:,k);
end

% Control hard constraints from 0 to N-1
c_ineq_u = nan(2*nu*N, 1);
for k = 1 : N
    c_ineq_u(2*nu*(k-1)+1 : 2*nu*(k-1)+2*nu) = [u_min-u_seq(:,k);
                                                u_seq(:,k)-u_max];
end

% Terminal constraint at N:
c_eq_xN = x_seq(:,N+1) - xt;

% Constraints on slack variables
c_ineq_s = nan(2*nx*(N+1),1);
for k = 1 : N+1
    c_ineq_s(2*nx*(k-1)+1 : 2*nx*(k-1)+2*nx) = -s_seq(:,k);
end

% % Concatenate constraints
% Inequality constraints
c_ineq = [c_ineq_x; % boundary constraint on state(from 0 to N)
          c_ineq_u; % boundary constraint on control(from 0 to N-1)
          c_ineq_s]; % non-negative constraint on slacks
          
% Equality constraint
c_eq = [c_eq_sc; % Continuity constraints (Multiple shooting method)
        c_eq_xN]; % terminal state == target setpoint

end

%%%%------------------------------%%%%------------------------------%%%%

%% Function helper: Cost function for DOA - 2 

function JN = getRegulationCostDOA(xi_seq, N)
% Dimensions
nx = 2; nu = 1;

% Regenerate slack variables from the decision variables
s_seq = xi_seq(nx*(N+1)+nu*N+1:end);

% Compute cost
JN = s_seq'*s_seq;

end

%%%%------------------------------%%%%------------------------------%%%%