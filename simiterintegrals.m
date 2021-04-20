function I = simiterintegrals(W, h, eps, varargin)
%SIMITERINTEGRALS Simulates twice iterated stochastic integrals.
%   I = SIMITERINTEGRALS(W,H,EPS) simulates the twice iterated stochastic
%   integral w.r.t the increment of a Wiener process W over stepsize H.
%   The algorithm guarantees that the error is less than or equal to EPS.
%   For SDEs the Maximum-L2 error estimate is used, whereas for SPDEs 
%   the Frobenius-L2 error estimate is used.
%
%   If the integrals are to be used in a SDE solver scheme then EPS should
%   typically be chosen as H^(p+0.5) where p is the order of the SDE scheme.
%   
%   I = SIMITERINTEGRALS(W,H,EPS,PARAM1,VAL1,PARAM2,VAL2,...) performs
%   the simulation with specified values of optional parameters.
%   The available parameters are
%
%   'ItoCorrection', whether to apply Ito correction
%
%       When 'ItoCorrection' is true the integrals are interpreted in the
%       Ito sense, otherwise in the sense of Stratonovich.
%       The default value is true.
%
%   'Algorithm', the algorithm to use for the simulation
%
%       Possible values are 'Fourier', 'Milstein', 'Wiktorsson' and 'MR'. 
%       The default value is 'MR'.
%
%   'MemUse', a scaling parameter for the memory usage
%
%       Possible values are in the interval [0,1]. This parameter scales
%       between the minimum required memory usage ('MemUse'=0) and the
%       maximum memory usage ('MemUse'=1). This is not a linear scale.
%       Typically a higher memory usage leads to faster execution times.
%       The default value is 1.
%
%   'q_12', eigenvalues of covariance of a Q-Wiener process with Q=diag(q_12)^2
%   
%       When 'q_12' is supplied, W is assumed to be a finite-dimensional
%       approximation of an infinite Q-Wiener process with covariance
%       matrix Q=diag(q_12)^2. 
%       Default value is a vector of ones.
%       
%

% dimension of the Wiener process
m = length(W);

% check input arguments
narginchk(3,inf);
ip = inputParser;
addRequired(ip,'WienerIncrement',@(x) isnumeric(x) && iscolumn(x));
addRequired(ip,'StepSize',@(x) (x>0) && isnumeric(x) && isscalar(x));
addRequired(ip,'Error',@(x) (x>0) && isnumeric(x) && isscalar(x));
addParameter(ip,'ItoCorrection',true,@(x) isscalar(x) && ...
    (islogical(x) || (isnumeric(x) && (x == 0 || x == 1))));
addParameter(ip,'Algorithm',"mr",@(x) isstring(x) || ischar(x));
addParameter(ip,'MemUse',1,@(x) isscalar(x) && isnumeric(x) && ...
    0<=x && x<=1);
addParameter(ip,'q_12',ones(m,1),@(x) isnumeric(x) && ...
    length(x)==m && iscolumn(x));
parse(ip,W,h,eps,varargin{:});

% calculate the number of terms needed for the given precision
% and the Levy area using the chosen algorithm
if all(ip.Results.q_12 == 1)
    % standard Wiener process
    switch lower(ip.Results.Algorithm)
        case "fourier"
            n = ceil( 1.5*(h/(pi*eps))^2 );
            I = h*levyarea_fourier(W/sqrt(h),n,ip.Results.MemUse);
        case "milstein"
            n = ceil( 0.5*(h/(pi*eps))^2 );
            I = h*levyarea_milstein(W/sqrt(h),n,ip.Results.MemUse);
        case "wiktorsson"
            n = ceil( sqrt(5*m/12) * h/(pi*eps) );
            I = h*levyarea_wik(W/sqrt(h),n,ip.Results.MemUse);
        case "mr"
            n = ceil( sqrt(m/12) * h/(pi*eps) );
            I = h*levyarea_mr(W/sqrt(h),n,ip.Results.MemUse);
        otherwise
            error("Unknown algorithm: " + ip.Results.Algorithm + ...
                ". Possible choices are: Fourier, Milstein, Wiktorsson, MR.");
    end
else
    % Q-Wiener process with covariance matrix Q=diag(q_12)^2
    sqcov = sqrt(h)*ip.Results.q_12;
    coeff = sqrt(sum(ip.Results.q_12.^2)^2-sum(ip.Results.q_12.^4));
    switch lower(ip.Results.Algorithm)
        case "fourier"
            n = ceil( 1.5*(coeff*h/(pi*eps))^2 );
            I = diag(sqcov)*levyarea_fourier(W./sqcov,n,ip.Results.MemUse)*diag(sqcov);
        case "milstein"
            n = ceil( 0.5*(coeff*h/(pi*eps))^2 );
            I = diag(sqcov)*levyarea_fourier(W./sqcov,n,ip.Results.MemUse)*diag(sqcov);
        case "wiktorsson"
            n = ceil( sqrt(5*m/12) * coeff*h/(pi*eps) );
            I = diag(sqcov)*levyarea_wik(W./sqcov,n,ip.Results.MemUse)*diag(sqcov);
        case "mr"
            n = ceil( sqrt(m/12) * coeff*h/(pi*eps) );
            I = diag(sqcov)*levyarea_mr(W./sqcov,n,ip.Results.MemUse)*diag(sqcov);
        otherwise
            error("Unknown algorithm for Q-Wiener processes: " + ...
                ip.Results.Algorithm + ...
                ". Possible choices are: Fourier, Milstein, Wiktorsson, MR.");
    end
end

% calculate iterated integrals
I = I + W.*W'./2;

% optionally add Ito correction
if ip.Results.ItoCorrection
    I = I - diag(ip.Results.q_12).^2.*h./2;
end

end % simiterintegrals


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