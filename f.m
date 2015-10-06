function xdot = f(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file                Author: Sean Phillips
%
% Project: Simulation of four coupled linear ocillators over a completely 
% connected network. 
%
% Description: Flow map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global N A B n
n = length(A);
% differential equations
xdot = kron(eye(N),A)*x(1:n*N) + kron(eye(N),B)*x(n*N+1: N*n+N);

etadot = zeros(N,1);
taudot = -1;
xdot = [xdot;etadot;taudot];

end