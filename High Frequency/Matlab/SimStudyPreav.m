%--------------------------------------------------------------------------
%  This code computes the bias and MSE of the pre-averaging estimator 
% with noisy high-frequency data and n=23400 returns, i.e. 1s sampling.
%--------------------------------------------------------------------------

rng('default'); %uses default random number generator. 
%Results (bias, MSE) depend slighlty on the random number generator
%Setting parameters
nsim = 10000;  %number of simulations
n = 23400;  % i.e. 1 sec sampling in 6.5 trading hours
theta = 1;  %tuning parameter of pre-averaging,
[sigma0,kappa,xi,rho] = deal(0.04/250,5/250,0.50/250,-0.50);
gamma = 0.5; % noise to signal ratio
%Simulation loop starts and P denotes the pre-averaging estimator
[P, IV] = deal(NaN(nsim,1));
for i = 1 : nsim
    fprintf('Sim no. %5d of %5d...\n',i,nsim); tic;
    [X,sigma] = f_SVHeston(1, n,sigma0,kappa,xi,rho); %Simulate Heston 
    IV(i) = mean(sigma(1:end-1).^2); %IV is approximated by a sum
     %Adding noise to prices
      omega = gamma*sqrt(IV(i)/n); %omega^2 is variance of the noise
      noise = omega*randn(n+1,1); %add iid gaussian noise
      Y =X+noise; %Y is the observed price
      % pre-averaging
      K= round(theta*sqrt(n)); %block size k_n
      P(i) = f_preav(Y,K);
      fprintf('Time elapsed: %5.2f...\n',toc);
end
bias=mean(P-IV);
Rbias=mean(P./IV-1);
MSE=mean((P-IV).^2);