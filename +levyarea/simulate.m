function I = simulate(W, h, eps, ito, alg, mem_use, q_12, err_norm)

% check input arguments
narginchk(8,8);

% dimension of the Wiener process
m = length(W);

% interpret error norm and set the scaling factor for Q-Wiener processes
% (Frobenius-L2 error < coeff * Maximum-L2 error)
switch lower(err_norm)
    case "maxl2"
        % maximum(q_12[i]*q_12[j] for i=1:m for j=1:i-1)
        Q = q_12*q_12';
        norm_coeff = max(Q(~eye(m)));
    case "frobeniusl2"
        norm_coeff = sqrt(sum(q_12.^2)^2-sum(q_12.^4));
    otherwise
        error("Unknown error norm: " + err_norm + ...
        ". Possible choices are: MaxL2, FrobeniusL2.");
end

% calculate the number of terms needed for the given precision
% and the Levy area using the chosen algorithm
sqcov = sqrt(h)*q_12;
switch lower(alg)
    case "fourier"
        n = ceil( 1.5*(norm_coeff*h/(pi*eps))^2 );
        I = diag(sqcov)*levyarea_fourier(W./sqcov,n,mem_use)*diag(sqcov);
    case "milstein"
        n = ceil( 0.5*(norm_coeff*h/(pi*eps))^2 );
        I = diag(sqcov)*levyarea_milstein(W./sqcov,n,mem_use)*diag(sqcov);
    case "wiktorsson"
        n = ceil( sqrt(5*m/12) * norm_coeff*h/(pi*eps) );
        I = diag(sqcov)*levyarea_wik(W./sqcov,n,mem_use)*diag(sqcov);
    case "mr"
        n = ceil( sqrt(m/12) * norm_coeff*h/(pi*eps) );
        I = diag(sqcov)*levyarea_mr(W./sqcov,n,mem_use)*diag(sqcov);
    otherwise
        error("Unknown algorithm for Q-Wiener processes: " + ...
            ip.Results.Algorithm + ...
            ". Possible choices are: Fourier, Milstein, Wiktorsson, MR.");
end

% calculate iterated integrals
I = I + W.*W'./2;

% optionally add Ito correction
if ito
    I = I - diag(q_12).^2.*h./2;
end

end % simulate


function I = levyarea_fourier(W, n, mem_use)
    m = length(W);
    
    if mem_use == 1
        alpha = randn(n,m);
        beta = (randn(m,n) - sqrt(2).*W) ./ (1:n);
        S = beta * alpha;
    elseif ~(0<=mem_use && mem_use<=1)
        error("Parameter 'mem_use' must be in [0,1].");
    else
        if mem_use == 0
            num_iter = n;
            n_small = 1;
        else
            num_iter = ceil(1/mem_use);
            n_small = ceil(n/num_iter);
        end
        S = zeros(m,m);
        for i = 1:num_iter
            alpha = randn(n_small,m);
            beta = (randn(m,n_small) - sqrt(2).*W) ./ ((i-1)*n_small+1:i*n_small);
            S = S + beta * alpha;
        end
    end
    
    I = (S-S')/(2*pi);
end % levyarea_fourier


function I = levyarea_milstein(W, n, mem_use)
    m = length(W);
    
    if mem_use == 1
        alpha = randn(n,m);
        beta = (randn(m,n) - sqrt(2).*W) ./ (1:n);
        S = beta * alpha;
        last_index = n;
    elseif ~(0<=mem_use && mem_use<=1)
        error("Parameter 'mem_use' must be in [0,1].");
    else
        if mem_use == 0
            num_iter = n;
            n_small = 1;
        else
            num_iter = ceil(1/mem_use);
            n_small = ceil(n/num_iter);
        end
        last_index = num_iter*n_small;
        S = zeros(m,m);
        for i = 1:num_iter
            alpha = randn(n_small,m);
            beta = (randn(m,n_small) - sqrt(2).*W) ./ ((i-1)*n_small+1:i*n_small);
            S = S + beta * alpha;
        end
    end
    
    M = randn(1,m);
    S = S + sqrt(2*psi(1,last_index+1)) * W.*M;
    I = (S-S')/(2*pi);
end % levyarea_milstein


function G = levyarea_wik(W, n, mem_use)
    m = length(W);
    
    if mem_use == 1
        alpha = randn(n,m);
        beta = (randn(m,n) - sqrt(2).*W) ./ (1:n);
        S = beta * alpha;
        last_index = n;
    elseif ~(0<=mem_use && mem_use<=1)
        error("Parameter 'mem_use' must be in [0,1].");
    else
        if mem_use == 0
            num_iter = n;
            n_small = 1;
        else
            num_iter = ceil(1/mem_use);
            n_small = ceil(n/num_iter);
        end
        last_index = num_iter*n_small;
        S = zeros(m,m);
        for i = 1:num_iter
            alpha = randn(n_small,m);
            beta = (randn(m,n_small) - sqrt(2).*W) ./ ((i-1)*n_small+1:i*n_small);
            S = S + beta * alpha;
        end
    end
    
    G = zeros(m);
    G(tril(true(m),-1)) = sqrt(2*psi(1,last_index+1)) .* randn(m*(m-1)/2,1);
    S = S + ((G-G')*W).*W'./(1+sqrt(1+norm(W)^2)) + G;
    G = (S-S')/(2*pi);
end % levyarea_wik


function G = levyarea_mr(W, n, mem_use)
    m = length(W);
    
    if mem_use == 1
        alpha = randn(n,m);
        beta = (randn(m,n) - sqrt(2).*W) ./ (1:n);
        S = beta * alpha;
        last_index = n;
    elseif ~(0<=mem_use && mem_use<=1)
        error("Parameter 'mem_use' must be in [0,1].");
    else
        if mem_use == 0
            num_iter = n;
            n_small = 1;
        else
            num_iter = ceil(1/mem_use);
            n_small = ceil(n/num_iter);
        end
        last_index = num_iter*n_small;
        S = zeros(m,m);
        for i = 1:num_iter
            alpha = randn(n_small,m);
            beta = (randn(m,n_small) - sqrt(2).*W) ./ ((i-1)*n_small+1:i*n_small);
            S = S + beta * alpha;
        end
    end
    
    M = randn(1,m);
    G = zeros(m);
    G(tril(true(m),-1)) = randn(m*(m-1)/2,1);
    S = S + sqrt(2*psi(1,last_index+1)) * (W.*M + G);
    G = (S-S')/(2*pi);
end % levyarea_mr