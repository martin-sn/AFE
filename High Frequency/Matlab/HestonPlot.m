%--------------------------------------------------------------------------
%   This code simulate Heston model and plots spot volatility, IV, RV
%   and confidence intervals for 10 days with 2340 returns in each day.
%                                                           
%--------------------------------------------------------------------------

% Setting parmeters
n = 2340; dt = 1/n; T=10; 
[sigma0,kappa,xi,rho] = deal(0.04/250,5/250,0.50/250,-0.50);
%Simulate X and sigma_t 
%(different random number generators may yield different X and sigma)
[X,sigma] = f_SVHeston(T, n,sigma0,kappa,xi,rho);
r=reshape(diff(X),n, T); % returns divided to T different days
%Realized volatility computed on different days
RV=sum(r.^2);
s=reshape(sigma(1:end-1), n, T);
IV = mean(s.^2); % IV is an intergal, which is approximated by a sum.


%Plotting sigma_t^2, RV and IV
a=100*250; %we multiply everything with number a for better visualization
figure; 
plot(0:dt:T,a*sigma.^2,'-b',1:1:T,a*RV,'kx',1:1:T, a*IV,'ro');
legend('Spot Vol', 'RV', 'IV');
xlabel('time');
ylabel('volatility');
hold off;
%Plotting RV with confidence intervals and IV for 10 days
figure;
RQ=n*sum(r.^4); %we need realized quarticity for CIs
CI=1.96*sqrt(2*RQ/3)/sqrt(n);
id=1:1:T;
plot(id, a*RV,'kx', id, a*IV,'ro', id, a*(RV-CI), 'kd', id, a*(RV+CI), 'kd');
legend('RV', 'IV', 'CI');
xlabel('time');
ylabel('volatility');
  