function [W, R_curr, P, rank_diff_all, R_all, P_all] = SCA_algorithm(para, h, rho_c, rho_r, reduction_factor)
%The successive convex approximization algorithm for NOMA-DFRC
%  [W,R_curr,P] = SCA_algorithm(para, h, rho)
%Inputs:
%   para: system paramaters
%   h: channel gains
%   rho: regularization parameter
%Outputs:
%   W: optimal beamformer
%   R_curr: optimal rate
%   P: optimal probing power
%Date: 14/10/2021
%Author: Zhaolin Wang


epsilon_inner = 1.001;
epsilon_outter = 1e-4;
eta = 1e5;

H = zeros(para.N, para.N, para.K);
for k = 1:para.K
    h_k = h(:,k);
    H(:,:,k) = h_k * h_k';
end

[W_n,cvx_status] = initial_point(para, H);
if strcmp(cvx_status, 'Infeasible')
    W = 0; R_curr = 0; P = 0;
    return;
end

%% To demostrate convergence
[P] = probing_power(para, sum(W_n,3));
[R_curr] = rate_NOMA(para, W_n, H);
rank_diff = 0;
for k = 1:para.K
    s = svd(W_n(:,:,k));
    rank_diff = rank_diff + sum(s) - max(s);
end

R_all = [sum(min(R_curr,[],2))];
P_all = [sum(P)];
rank_diff_all = [rank_diff];

%% SCA algorithm
rank_diff = 100;
f_pre = -1;
while rank_diff > epsilon_outter
    f_frac = 100;
	while f_frac > epsilon_inner

        [W,f,P,cvx_status] = SCA_step(para, W_n, H, rho_c, rho_r, eta);
        
        if strcmp(cvx_status, 'Failed') || strcmp(cvx_status, 'Infeasible')
            W = 0; R_curr = 0; P = 0;
            return;
        else
            [R_curr] = rate_NOMA(para, W, H);
            disp(['sum rate -- ', num2str(sum(min(R_curr,[],2))),' / sum probing power-- ', num2str(sum(P)) ' / objective value', num2str(f)]);
            f_frac = min(f, f_pre) / max(f, f_pre);
            f_pre = f;
            W_n = W;
        end
    end
    
    rank_diff = 0;
    for k = 1:para.K
        s = svd(W_n(:,:,k));
        rank_diff = rank_diff + sum(s) - max(s);
    end
    eta = reduction_factor * eta;
    disp(['penalty term -- ', num2str(rank_diff)]);
    
    % To demostrate convergence
    rank_diff_all = [rank_diff_all, rank_diff];
    R_all = [R_all, sum(min(R_curr,[],2))];
    P_all = [P_all, sum(P)];

end
[R_curr] = rate_NOMA(para, W, H);
end


function [W,f,P,cvx_status] = SCA_step(para, W_n, H, rho_c, rho_r, eta)
    cvx_begin quiet
    
        variable W(para.N,para.N,para.K) complex
        variable r(para.K, 1)
        Rw = sum(W,3);
        [P] = probing_power(para, Rw);
    
        [sca] = SCA_SINR(para, H, W, W_n, r);
    
        [corr] = cross_correlation(para, Rw);
        
        % constraints
        for k = 1:para.K
            W(:,:,k) == hermitian_semidefinite(para.N);
            r(k) >= para.R_min;
        end
    
        for i = 1:length(sca)
            sca(i) >= 0;
        end
        M = length(para.theta);
        for k = 1:M
            for p = k+1:M
                P(k) - P(p) <= 10;
                P(k) - P(p) >= -10;
            end
        end
    
        diag(Rw) == para.Pt / para.N;
        corr <= para.correlation_min;
        
        % objective function
        [f] = obj_func(W, W_n, r, P, para.K, rho_c, rho_r, eta);
        minimize (f);
    cvx_end
end



%% objective function
function [f] = obj_func(W, W_n, r, P, K, rho_c, rho_r, eta)
    f = -rho_c*sum(r) - rho_r*sum(P) + 1/eta*rank_penalty(W, W_n, K);
end

%% rank-1 penalty term
function [p] = rank_penalty(W, W_n, K)
    p = 0;
    for k = 1:K
        W_k = W(:,:,k);
        W_n_k = W_n(:,:,k);
        
        [~,S,V] = svd(W_n_k);
        v_max = V(:,1);
        spectral_taylor = -S(1,1) - real(trace(v_max * v_max' * (W_k - W_n_k)));
        
        p = p + norm_nuc(W_k) + spectral_taylor;
    end
end

%% initialize
function [W,cvx_status] = initial_point(para, H)

    cvx_begin quiet
        
        variable W(para.N,para.N,para.K) complex

        Rw = sum(W,3);
        [P] = probing_power(para, Rw);
        [ineqn] = initial_rate(para,H, W);
        [corr] = cross_correlation(para, Rw);
        for k = 1:para.K
            W(:,:,k) == hermitian_semidefinite(para.N);
        end
        
        for i = 1:length(ineqn)
            ineqn(i) <= 0;
        end
        
        M = length(para.theta);
        for k = 1:M
            for p = k+1:M
                P(k) - P(p) <= 10;
                P(k) - P(p) >= -10;
            end
        end

        diag(Rw) == para.Pt / para.N;
        corr <= para.correlation_min;
    cvx_end   
end
