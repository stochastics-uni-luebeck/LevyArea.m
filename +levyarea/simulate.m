function I = simulate(W, h, err, ito, alg, q_12, err_norm)
%SIMULATE The core routine to simulate iterated integrals.
% See also ITERATED_INTEGRALS, OPTIMAL_ALGORITHM.

% check number of input arguments
narginchk(7,7);

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
        n = ceil( 1.5*(norm_coeff*h/(pi*err))^2 );
        I = diag(sqcov)*levyarea_fourier(W./sqcov,n)*diag(sqcov);
    case "milstein"
        n = ceil( 0.5*(norm_coeff*h/(pi*err))^2 );
        I = diag(sqcov)*levyarea_milstein(W./sqcov,n)*diag(sqcov);
    case "wiktorsson"
        n = ceil( sqrt(5*m/12) * norm_coeff*h/(pi*err) );
        I = diag(sqcov)*levyarea_wik(W./sqcov,n)*diag(sqcov);
    case "mr"
        n = ceil( sqrt(m/12) * norm_coeff*h/(pi*err) );
        I = diag(sqcov)*levyarea_mr(W./sqcov,n)*diag(sqcov);
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


function I = levyarea_fourier(W, n)
    m = length(W);
    alpha = randn(n,m);
    beta = (randn(m,n) - sqrt(2).*W) ./ (1:n);
    S = beta * alpha;
    I = (S-S')/(2*pi);
end % levyarea_fourier


function I = levyarea_milstein(W, n)
    m = length(W);
    alpha = randn(n,m);
    beta = (randn(m,n) - sqrt(2).*W) ./ (1:n);
    S = beta * alpha;
    M = randn(1,m);
    S = S + sqrt(2*psi(1,n+1)) * W.*M;
    I = (S-S')/(2*pi);
end % levyarea_milstein


function G = levyarea_wik(W, n)
    m = length(W);
    alpha = randn(n,m);
    beta = (randn(m,n) - sqrt(2).*W) ./ (1:n);
    S = beta * alpha;
    G = zeros(m);
    G(tril(true(m),-1)) = sqrt(2*psi(1,n+1)) .* randn(m*(m-1)/2,1);
    S = S + ((G-G')*W).*W'./(1+sqrt(1+norm(W)^2)) + G;
    G = (S-S')/(2*pi);
end % levyarea_wik


function G = levyarea_mr(W, n)
    m = length(W);
    alpha = randn(n,m);
    beta = (randn(m,n) - sqrt(2).*W) ./ (1:n);
    S = beta * alpha;
    M = randn(1,m);
    G = zeros(m);
    G(tril(true(m),-1)) = randn(m*(m-1)/2,1);
    S = S + sqrt(2*psi(1,n+1)) * (W.*M + G);
    G = (S-S')/(2*pi);
end % levyarea_mr

