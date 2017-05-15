
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

%% Part 1
%%%%%%%%%%%%%%%%%%%%%    First MPC controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART I - First MPC controller...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup
n_states = 7;
n_input = 4;
N = 5; % TODO: Horizon length does not affect system output. Is this correct??
N = N+1; % account for the fact that mpc is defined for x0-xN, but matlab array indexing starts at 1.

% cost matrices
Q = diag([0.5,10,10,0.5,0.5,0.5,0.5])
R = 0.1*eye(n_input)


% yalmip variables
x = sdpvar(n_states, N);
u = sdpvar(n_input, N);

% state constraints
z_dot_max = 1;
alpha_beta_max = 10/360 * 2 * pi;
alpha_beta_dot_max = 15 / 360 * 2 * pi;
gamma_dot_max = 60 / 360 * 2 * pi;
u_min = 0-us;
u_max = 1-us;

x0 =  [-1 0.1745 -0.1745 0.8727 0 0 0]';

% Terminal Set

% Terminal set it positive invariant set for the system x(k+1) = (A +
% BF)x(k), where F is constructed as the infinite horizon control input and
% Pinf is the terminal constraint.

[Pinf,L,Finf] = dare(sys.A,sys.B,Q,R);

terminal_model = LTISystem('A',sys.A - sys.B * Finf);

mpt_mpc = MPCController(terminal_model,N-1);

fprintf('calculate invariant set\n')
terminal_set = mpt_mpc.model.invariantSet();

% deliverables Ax and bx
terminal_set.A
terminal_set.b

terminal_constraint = ismember(x(:,N),terminal_set);


% Constraints
fprintf('setting constraints \n')
constraints_mpc = [terminal_constraint];
for i = 2:N
    % set system constraints
    constraints_mpc = constraints_mpc + [x(:,i) == sys.A* x(:,i-1) + sys.B * u(:,i-1)];
    
    % set state constraints
    constraints_mpc =  constraints_mpc +  [-z_dot_max <= x(1,i) <= z_dot_max];
    constraints_mpc =  constraints_mpc + [-alpha_beta_max <= x(2:3,i) <= alpha_beta_max];
    constraints_mpc = constraints_mpc + [-alpha_beta_dot_max <= x(5:6,i) <= alpha_beta_dot_max];
    constraints_mpc = constraints_mpc + [-gamma_dot_max <= x(7,i) <= gamma_dot_max];

    % set input constraints
    constraints_mpc = constraints_mpc + [u_min <= u(1:4,i) <= u_max];
end
constraints_mpc = constraints_mpc + [u_min <= u(1:4,1) <= u_max]; % constrain u(1) alias u0 as well.

% Objective / Cost Function
fprintf('setting objective function\n')
objective_mpc = 0;
for i = 1:N-1
    objective_mpc = objective_mpc + x(:,i)' * Q * x(:,i) + u(:,i)' * R * u(:,i);
end

terminal_cost = x(:,N)' * Pinf * x(:,N);
objective_mpc = objective_mpc + terminal_cost;

% set first column to x0
parameters = x(:,1);
wanted = u(:,1);

% Setup MPC
fprintf('set up mpc problem\n')
mpc_simple = optimizer(constraints_mpc, objective_mpc, [], parameters, wanted);

% Simulation
fprintf('simulate system\n')
simQuad( sys, mpc_simple, 0, x0, 20);

%% Part 2
%%%%%%%%%%%%%%%%%%%%%  Reference Tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART II - Reference tracking...\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup
n_states = 7;
n_input = 4;
N =5; 
N = N+1; % account for the fact that mpc is defined for x0-xN, but matlab array indexing starts at 1.

% cost matrices
Q = diag([1,50,50,1,0.1,0.1,0.1])
R = 0.1*eye(n_input)
P = eye(n_states)
Rs = eye(n_input)

% yalmip variables
delta_x = sdpvar(n_states, N,'full');
delta_u = sdpvar(n_input, N,'full');

r = sdpvar(4,1,'full');
xr = sdpvar(n_states,1,'full')
ur = sdpvar(n_input,1,'full')

xk = sdpvar(n_states,1,'full')
uk = sdpvar(n_input,1,'full')

C = [eye(4) zeros(4,3)];

% state constraints
z_dot_max = 1;
alpha_beta_max = degtorad(10);
alpha_beta_dot_max = degtorad(15) ;
gamma_dot_max = degtorad(60);
u_min = 0-us
u_max = 1-us

x0 =  [-1 0.1745 -0.1745 0.8727 0 0 0]';
r1 = [1 0.1745 -0.1745 1.7453]'; 

% compute target state ts: (xr, ur)
xrur = inv([eye(n_states)-sys.A -sys.B; C zeros(4,n_input)])*[zeros(n_states,1); r1]
%xr = [r1; 0; 0; 0];
%ur = xrur(n_states + 1 : end);

% Constraints
fprintf('setting constraints \n')

constraints_mpc = [];

%target constraints
constraints_mpc = constraints_mpc + [C*xr == r];
constraints_mpc = constraints_mpc + [sys.A*xr+sys.B*ur == xr];

%delta shifting
constraints_mpc = constraints_mpc + [delta_x(:,1) == xk-xr];
constraints_mpc = constraints_mpc + [delta_u(:,1) == uk-ur];
 
for i = 2:N
    %system constraints
    constraints_mpc = constraints_mpc + [delta_x(:,i) == sys.A*delta_x(:,i-1)+sys.B*delta_u(:,i-1)];
    
    %state constraints
    constraints_mpc =  constraints_mpc +  [-z_dot_max-xr(1) <= delta_x(1,i) <= z_dot_max-xr(1)];
    constraints_mpc =  constraints_mpc + [-alpha_beta_max-xr(2:3)  <= delta_x(2:3,i) <= alpha_beta_max-xr(2:3)];
    constraints_mpc = constraints_mpc + [-alpha_beta_dot_max-xr(5:6)  <= delta_x(5:6,i) <= alpha_beta_dot_max-xr(5:6)];
    constraints_mpc = constraints_mpc + [-gamma_dot_max-xr(7)  <= delta_x(7,i) <= gamma_dot_max-xr(7)];

    %input constraints
    
    constraints_mpc = constraints_mpc + [u_min-ur <= delta_u(1:4,i) <= u_max-ur ];
end
constraints_mpc = constraints_mpc + [u_min-ur <= delta_u(1:4,1) <= u_max-ur ];

%Objective / Cost Function
fprintf('setting objective function\n')
objective_mpc = 0;
for i = 1:N
    objective_mpc = objective_mpc + delta_x(:,i)' * Q * delta_x(:,i) + delta_u(:,i)' * R * delta_u(:,i);
end

[Pinf,L,Finf] = dare(sys.A,sys.B,Q,R);

terminal_cost = delta_x(:,N)' * Pinf * delta_x(:,N);
objective_mpc = objective_mpc + terminal_cost;

% Setup MPC
fprintf('set up mpc problem\n')

mpc_simple = optimizer(constraints_mpc, objective_mpc, [], [xk; r], uk);


% Simulation constant r
fprintf('simulate system with constant r\n')
simQuad( sys, mpc_simple, 0, zeros(7,1), 10, r1);

fprintf('simulate system with varying r r\n')
n_k = 150
k = 1:n_k;
rt = [repmat(1,1,n_k); 0.1745*sin(sys.Ts*k); -0.1745*sin(sys.Ts*k); repmat(pi/2,1,n_k) ];
simQuad( sys, mpc_simple, 0, zeros(7,1), 20, rt);


%% %%%%%%%%%%%%%  First simulation of the nonlinear model %%%%%%%%%%%%%%%%%
fprintf('PART III - First simulation of the nonlinear model...\n')
close all; 
innerController = mpc_simple; 
sim('simulation1.mdl'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%  Offset free MPC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('PART IV - Offset free MPC...\n')

% define L
Lx = eye(7);
Ld = eye(7); 
L = [Lx; Ld]; 

% Setup
n_states = 7;
n_dist = 7; 
n_input = 4;
N =5; 
N = N+1; % account for the fact that mpc is defined for x0-xN, but matlab array indexing starts at 1.

% cost matrices
Q = diag([1,50,50,1,0.1,0.1,0.1])
R = 0.1*eye(n_input)
P = eye(n_states)
Rs = eye(n_input)

% yalmip variables
delta_x = sdpvar(n_states, N,'full');
delta_u = sdpvar(n_input, N,'full');

r = sdpvar(4,1,'full');
xr = sdpvar(n_states,1,'full')
ur = sdpvar(n_input,1,'full')

xk = sdpvar(n_states,1,'full')
uk = sdpvar(n_input,1,'full')

C = [eye(4) zeros(4,3)];

% state constraints
z_dot_max = 1;
alpha_beta_max = degtorad(10);
alpha_beta_dot_max = degtorad(15) ;
gamma_dot_max = degtorad(60);
u_min = 0-us
u_max = 1-us

x0 =  [-1 0.1745 -0.1745 0.8727 0 0 0]';
r1 = [1 0.1745 -0.1745 1.7453]'; 

% compute target state ts: (xr, ur)
xrur = inv([eye(n_states)-sys.A -sys.B; C zeros(4,n_input)])*[zeros(n_states,1); r1]
%xr = [r1; 0; 0; 0];
%ur = xrur(n_states + 1 : end);

% Constraints
fprintf('setting constraints \n')

constraints_mpc = [];

%target constraints
constraints_mpc = constraints_mpc + [C*xr == r];
constraints_mpc = constraints_mpc + [sys.A*xr+sys.B*ur == xr];

%delta shifting
constraints_mpc = constraints_mpc + [delta_x(:,1) == xk-xr];
constraints_mpc = constraints_mpc + [delta_u(:,1) == uk-ur];
 
for i = 2:N
    %system constraints
    constraints_mpc = constraints_mpc + [delta_x(:,i) == sys.A*delta_x(:,i-1)+sys.B*delta_u(:,i-1)];
    
    %state constraints
    constraints_mpc =  constraints_mpc +  [-z_dot_max-xr(1) <= delta_x(1,i) <= z_dot_max-xr(1)];
    constraints_mpc =  constraints_mpc + [-alpha_beta_max-xr(2:3)  <= delta_x(2:3,i) <= alpha_beta_max-xr(2:3)];
    constraints_mpc = constraints_mpc + [-alpha_beta_dot_max-xr(5:6)  <= delta_x(5:6,i) <= alpha_beta_dot_max-xr(5:6)];
    constraints_mpc = constraints_mpc + [-gamma_dot_max-xr(7)  <= delta_x(7,i) <= gamma_dot_max-xr(7)];

    %input constraints
    
    constraints_mpc = constraints_mpc + [u_min-ur <= delta_u(1:4,i) <= u_max-ur ];
end
constraints_mpc = constraints_mpc + [u_min-ur <= delta_u(1:4,1) <= u_max-ur ];

%Objective / Cost Function
fprintf('setting objective function\n')
objective_mpc = 0;
for i = 1:N
    objective_mpc = objective_mpc + delta_x(:,i)' * Q * delta_x(:,i) + delta_u(:,i)' * R * delta_u(:,i);
end

[Pinf,L,Finf] = dare(sys.A,sys.B,Q,R);

terminal_cost = delta_x(:,N)' * Pinf * delta_x(:,N);
objective_mpc = objective_mpc + terminal_cost;

% Setup MPC
fprintf('set up mpc problem\n')

mpc_simple = optimizer(constraints_mpc, objective_mpc, [], [xk; r], uk);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%  Simulation of the nonlinear model %%%%%%%%%%%%%%%%%%%%
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
