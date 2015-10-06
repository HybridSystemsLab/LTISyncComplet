function xplus = g(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file               Author: Sean Phillips
%
% Project: Simulation of four coupled linear ocillators over a completely 
% connected network. 
%
% Description: Jump map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global K N G v n B

xg = x(1:N*n);
p = size(B);
etaplus = NaN(N*p(2),1);

for i = 1:N
    e = kron(G(:,i),eye(n))'*(kron(ones(N,1),eye(n))*x((i-1)*n+1:i*n) - xg);
    etaplus(i) = K*1/N*e;
end

tauplus = v(2) - (v(2) - v(1))*rand(1);
xplus = [x(1:n*N);etaplus;tauplus];

end