% initial conditions
clear all
close all

% Used Only for scalar systems

global K v N G A B n
A = [0 1; -1 0];
B = [0; 1];
n = length(A);
N = 4;
G = [0 1 1 1; 1 0 1 1; 1 1 0 1; 1 1 1 0];
% G = [0 1 1 0; 1 0 0 1; 1 0 0 1; 1 0 0 0]; 
K = [-.5, -.7];
v = [0.3 .65];

% Initiallizing Parameters and Intial Conditions
n = length(A); %Size of system matrix
N = length(G); %Number of agents
% Initial conditions for each agent 
x10 = [-1 0]; % Agent 1
eta10 = 1.76;
tau10 = 0.18;
x20 = [1 0]; % Agent 2
eta20 = 1.98;
tau20 = 0.29;
x30 = [0.5 0]; % Agent 3
eta30 = 0; 
tau30 = 0.15;
x40 = [2 0]; % Agent 4
eta40 = 1.73;
tau40 = 0.14;

x0 = [x10 x20 x30 x40];
eta0 = [eta10 eta20 eta30 eta40];
tau0 = [tau10];
X0 = [x0 eta0 tau0]';
% simulation horizon
TSPAN=[0 60];
JSPAN = [0 1000];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.1);

% simulate hybrid system using hybridsolver.m
[t, y, j] = hybridsolver( @f,@g,@C,@D,X0,TSPAN,JSPAN,rule,options,1);

% Using parameters P and sigma given by 
sigma = 0.1;
P1 = [23.6327001603891,5.10984623811966,3.29719921773053;5.10984623811966,13.6056403999659,1.84052617306014;3.29719921773053,1.84052617306014,7.51232754363112];

% Define Af, Bf, Ag for oscilator system
Af = [A B; 0 0 0];
Bf = [zeros(size(A)), B; zeros(size(B')), 0];
Ag = [eye(size(A)), zeros(size(B)); K, K*zeros(size(B))];

% Evaluate Lyapunov function over the generated solutions 
Vi = NaN(1,length(t));
V = 0;
xg = y(:,1:n*N)';
etag = y(:,n*N+1:n*N+N)'; 
for i =1:N
    tau = y(:,N*n+N+1);
    
    for j = 1:length(t)
        err = 1/3*(4*xg((i-1)*n+1:(i*n),j) - xg(1:2,j) - xg(3:4,j) - xg(5:6,j) - xg(7:8,j));
        etaerr = 1/3*(4*etag(i,j) - etag(1,j) - etag(2,j) - etag(3,j) - etag(4,j));
        agent = [err',etaerr];
        Vi(j) = exp(sigma*tau(j))*agent*expm(Af'*tau(j))*P1*expm(Af*tau(j))*agent';
    end
    subplot(313)
    hold on
    V = V + Vi;
end

% Generate plots of data

x11 = y(:,1);
x12 = y(:,2);
x21 = y(:,3);
x22 = y(:,4);
x31 = y(:,5);
x32 = y(:,6);
x41 = y(:,7);
x42 = y(:,8);
figure(1)
set(1,'Position',[212 888 560 209])
subplot(1,2,2)
plot(t, V, 'Color', [.47 .67 .19]);
axis([0, 10, 0 230])
grid on
% axis([0,8,0,300])
subplot(1,2,1)
plot(x11,x12,'-.',x21,x22,':',x31,x32,x41,x42,'--')
hold on
plot(x11(1),x12(1),'*',...
    'MarkerEdgeColor',[0,.45,.74])
plot(x21(1),x22(1),'*',...
    'MarkerEdgeColor',[0.85,.33,.1])
plot(x31(1),x32(1),'*',...
    'MarkerEdgeColor',[0.93,.69,.13])
plot(x41(1),x42(1),'*',...
    'MarkerEdgeColor',[0.49,.18,.56])
grid on
legend('', '', '', '')
