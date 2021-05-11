function [X,sigma] = f_SVHeston(T,n,theta,kappa,xi,rho)
%--------------------------------------------------------------------------
%   Simulate Heston model
%     T:     num of days
%     n:     number of observations wuthin a day
%     theta: average variance
%     kappa: speed of mean reversion
%     xi: volatility-of-volatility
%     rho:   leverage correlation
%                                                           
%--------------------------------------------------------------------------
%Setting horizon and discretization
N = T*n + 1;
dt = 1/n;

%Initialization
v0=theta; %This always starts sigma_t^2 from theta
%shape = 2*theta*kappa/(xi^2);%Better option might be to start from the  
%scale = 1/(2*kappa/(xi^2)); %stationary distrubution of sigma_t^2 
%v0 = gamrnd(shape,scale); %which is gamma

% simulate v_t=sigma_t^2 which is a square root diffusion process
dW= randn(N,1)*sqrt(dt); % increments of Brownian motion W
v = ones(N,1)*v0;
dB= rho*dW + sqrt((1-rho^2))*randn(N,1)*sqrt(dt); % increments of B
for j = 2 : N
v(j)= v(j-1)+kappa*(theta - v(j-1))*dt+ xi*sqrt(v(j-1))*dB(j-1);
end
   
v=max(v,0); % to avoid negative values caused by the discretization.
sigma = sqrt(v);
X= cumsum(sigma.*dW);

end






