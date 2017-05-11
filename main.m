
clear all
close all
clc

% adding the subfolders to the path
addpath(genpath('functions'))
addpath(genpath('data'))

% loads:
%    hovering equilibrium (xs,us)
%    continuous time matrices Ac,Bc of the linearization
%    matrices sys.A, sys.B of the inner-loop discretized with sampling period sys.Ts
%    outerController optimizer instance for the outer controller
load('quadData.mat')
outerController = getOuterController(Ac);
disp('Data successfully loaded')

%%%%%%%%%%%%%%%% ADD YOUR CODE BELOW THIS LINE %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%    First MPC controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART I - First MPC controller...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup
n_states = 7;
n_input = 4;
N = 5;
N = N+1; % account for the fact that mpc is defined for x0-xN, but matlab array indexing starts at 1.


Q = diag([1,10,10,0.5,1,1,1]);
R = 0.1*eye(n_input);
P = eye(n_states);

%variables
x = sdpvar(n_states, N);
u = sdpvar(n_input, N);

% state constraints
z_max = 1;
alpha_beta_max = 10/360 * 2 * pi;
alpha_beta_dot_max = 15 / 360 * 2 * pi;
gamma_dot_max = 60 / 360 * 2 * pi;
u_min = 0-us;
u_max = 1-us;

x0 =  [-1 0.1745 -0.1745 0.8727 0 0 0]';

%% Terminal Set

% Terminal set it positive invariant set for the system x(k+1) = (A +
% BF)x(k), where F is constructed as the infinite horizon control input and
% Pinf is the terminal constraint.

[Pinf,L,Finf] = dare(sys.A,sys.B,Q,R);

terminal_model = LTISystem('A',sys.A - sys.B * Finf, 'Ts', sys.Ts);

mpt_mpc = MPCController(terminal_model,N-1);
fprintf('calculate invariant st\n')
terminal_set = mpt_mpc.model.invariantSet()

terminal_constraint = ismember(x(:,N),terminal_set);


%% Constraints
%set constraints. order does not matter.
fprintf('setting constraints \n')
constraints = [terminal_constraint];
for i = 2:N
    % set system constraints
    constraints = constraints + [x(:,i) == sys.A* x(:,i-1) + sys.B * u(:,i-1)];
    
    % set state constraints
    constraints =  constraints +  [-z_max <= x(2:3,i) <= z_max];
    constraints =  constraints + [-alpha_beta_max <= x(2:3,i) <= alpha_beta_max];
    constraints = constraints + [-alpha_beta_dot_max <= x(5:6,i) <= alpha_beta_dot_max];
    constraints = constraints + [-gamma_dot_max <= x(5:6,i) <= gamma_dot_max];

    % set input constraints
    constraints = constraints + [u_min <= u(1:4,i) <= u_max];
end
constraints = constraints + [u_min <= u(1:4,1) <= u_max]; % constrain u(1) alias u0 as well.

% define objective function
fprintf('setting objective function\n')
objective = 0;
for i = 1:N-1
    objective = objective + x(:,i)' * Q * x(:,i) + u(:,i)' * R * u(:,i);
end

terminal_cost = x(:,N)' * Pinf * x(:,N);
objective = objective + terminal_cost;

% set first column to x0
parameters = x(:,1);
wanted = u(:,1);

%% Do MPC
fprintf('set up mpc problem\n')
mpc_simple = optimizer(constraints, objective, [], parameters, wanted);

u_real = zeros(n_input,N);
x_real = zeros(n_states,N);

x_real(:,1) = x0;

fprintf('do mpc\n')

for i=1:N-1
    u_real(:,i)= mpc_simple(x_real(:,i));
    x_real(:,i+1) = sys.A * x_real(:,i) + sys.B * u_real(:,i);
end

%% simulation
simQuad( sys, mpc_simple, 0, x0, 20);

%%%%%%%%%%%%%%%%%%%%%  Reference Tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART II - Reference tracking...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%  First simulation of the nonlinear model %%%%%%%%%%%%%%%%%
fprintf('PART III - First simulation of the nonlinear model...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%  Offset free MPC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART IV - Offset free MPC...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%  Simulation of the nonlinear model %%%%%%%%%%%%%%%%%%%%
fprintf('PART V - simulation of the nonlinear model...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%  Slew Rate Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VI - Slew Rate Constraints...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%  Soft Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VII - Soft Constraints...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FORCES Pro %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART VIII - FORCES Pro...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
