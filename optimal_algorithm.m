function alg = optimal_algorithm(dim, h, eps, varargin)
%OPTIMAL_ALGORITHM Determines the optimal algorithm for SIMITERINTEGRALS.
%   ALG = OPTIMAL_ALGORITHM(DIM, H, EPS) determines the optimal algorithm
%   under the given parameters, i.e. the algorithm that needs to simulate
%   the fewest random numbers to achieve the desired precision.
%
%   This function accepts the optional parameters 'q_12' and 'ErrorNorm' as
%   specified for SIMITERINTEGRALS.

% check input arguments
narginchk(3,inf);
ip = inputParser;
addRequired(ip,'Dimension',@(x) (x>0) && isnumeric(x) && isscalar(x));
addRequired(ip,'StepSize',@(x) (x>0) && isnumeric(x) && isscalar(x));
addRequired(ip,'Error',@(x) (x>0) && isnumeric(x) && isscalar(x));
addParameter(ip,'q_12',ones(dim,1),@(x) isnumeric(x) && ...
    length(x)==dim && iscolumn(x));
addParameter(ip,'ErrorNorm',"auto",@(x) isstring(x) || ischar(x))
parse(ip,dim,h,eps,varargin{:});

% interpret error norm and set the scaling factor for Q-Wiener processes
% (Frobenius-L2 error < coeff * Maximum-L2 error)
switch lower(ip.Results.ErrorNorm)
    case "maxl2"
        % maximum(q_12[i]*q_12[j] for i=1:m for j=1:i-1)
        Q = ip.Results.q_12*ip.Results.q_12';
        norm_coeff = max(Q(~eye(m)));
    case "frobeniusl2"
        norm_coeff = sqrt(sum(ip.Results.q_12.^2)^2-sum(ip.Results.q_12.^4));
    case "auto"
        if all(ip.Results.q_12 == 1)
            % maximum(q_12[i]*q_12[j] for i=1:m for j=1:i-1)
            Q = ip.Results.q_12*ip.Results.q_12';
            norm_coeff = max(Q(~eye(m)));
        else
            norm_coeff = sqrt(sum(ip.Results.q_12.^2)^2-sum(ip.Results.q_12.^4));
        end
    otherwise
        error("Unknown error norm: " + ip.Results.ErrorNorm + ...
        ". Possible choices are: Auto, MaxL2, FrobeniusL2.");
end

% available Algorithms
algs = ["Fourier","Milstein","Wiktorsson","MR"];

% Fourier
n = ceil( 1.5*(norm_coeff*h/(pi*eps))^2 );
costs(1) = 2*n*dim;
% Milstein
n = ceil( 0.5*(norm_coeff*h/(pi*eps))^2 );
costs(2) = 2*n*dim+dim;
% Wiktorsson
n = ceil( sqrt(5*dim/12) * norm_coeff*h/(pi*eps) );
costs(3) = 2*n*dim+(dim^2-dim)/2;
% MR
n = ceil( sqrt(dim/12) * norm_coeff*h/(pi*eps) );
costs(4) = 2*n*dim+(dim^2-dim)/2+dim;

% find minimum
[~,idx] = min(costs);
alg = algs(idx);