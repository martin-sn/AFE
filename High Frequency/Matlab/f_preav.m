%--------------------------------------------------------------------------
%   Pre-averaging estimator
%       Y:       observed prices
%       K:       pre-averaging window, i.e. k_n
%       n:       sample size
%--------------------------------------------------------------------------
function P = f_preav(Y,K)
    kernel = @(x) min(x,1-x);
    r = diff(Y); % returns of observed prices 
    n =size(r,1);
   % we first find pre-averaged retunrs, i.e. \bar{Y}_{i/n}
	r_pa = 0;
    for i = 1 : K - 1
        r_pa = r_pa + kernel(i/K)*r(i:end+i-(K-1),:);
    end
    % psi1 and psi2 functions
        idx = (0:K)/K; g = kernel(idx); dg = diff(g);
        psi1 = sum(dg.^2)*K; 
        psi2 = sum(g.^2)/K;
    %Finding the final estimator
    P= sum(r_pa.^2);
    P=P*n/((n-K+2)*K*psi2);
    P=P-sum(r.^2)*psi1/(2*psi2*K^2);
end


