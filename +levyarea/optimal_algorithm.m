function alg = optimal_algorithm(dim, h, err, varargin)
%OPTIMAL_ALGORITHM Determines the optimal algorithm for iterated_integrals.
%   Alg = OPTIMAL_ALGORITHM(Dim,H,Err) determines the optimal algorithm
%   under the given parameters, i.e. the algorithm that needs to simulate
%   the fewest random numbers to achieve the desired precision.
%
%   Alg = OPTIMAL_ALGORITHM(___,Name,Value)
%   This function accepts the optional parameters 'QWiener' and 'ErrorNorm' as
%   specified for iterated_integrals.
%
%   See also ITERATED_INTEGRALS.

% check input arguments
narginchk(3,inf);
ip = inputParser;
addRequired(ip,'Dimension',@(x) (x>0) && isnumeric(x) && isscalar(x));
addRequired(ip,'StepSize',@(x) (x>0) && isnumeric(x) && isscalar(x));
addRequired(ip,'Error',@(x) (x>0) && isnumeric(x) && isscalar(x));
addParameter(ip,'QWiener',ones(dim,1),@(x) isnumeric(x) && ...
    length(x)==dim && iscolumn(x));
addParameter(ip,'ErrorNorm',"auto",@(x) isstring(x) || ischar(x))
parse(ip,dim,h,err,varargin{:});

if ip.Results.ErrorNorm == "auto"
    if all(ip.Results.QWiener == 1)
        err_norm = "maxl2";
    else
        err_norm = "frobeniusl2";
    end
else
    err_norm = ip.Results.ErrorNorm;
end

% interpret error norm and set the scaling factor for Q-Wiener processes
% (Frobenius-L2 error < coeff * Maximum-L2 error)
q_12 = ip.Results.QWiener;
switch lower(err_norm)
    case "maxl2"
        % maximum(q_12[i]*q_12[j] for i=1:m for j=1:i-1)
        Q = q_12*q_12';
        norm_coeff = max(Q(~eye(dim)));
    case "frobeniusl2"
        norm_coeff = sqrt(sum(q_12.^2)^2-sum(q_12.^4));
    otherwise
        error("Unknown error norm: " + err_norm + ...
        ". Possible choices are: Auto, MaxL2, FrobeniusL2.");
end

% available Algorithms
algs = ["Fourier","Milstein","Wiktorsson","MronRoe"];

% Fourier
n = ceil( 1.5*(norm_coeff*h/(pi*err))^2 );
costs(1) = 2*n*dim;
% Milstein
n = ceil( 0.5*(norm_coeff*h/(pi*err))^2 );
costs(2) = 2*n*dim+dim;
% Wiktorsson
n = ceil( sqrt(5*dim/12) * norm_coeff*h/(pi*err) );
costs(3) = 2*n*dim+(dim^2-dim)/2;
% Mrongowius-Roessler
n = ceil( sqrt(dim/12) * norm_coeff*h/(pi*err) );
costs(4) = 2*n*dim+(dim^2-dim)/2+dim;

% find minimum
[~,idx] = min(costs);
alg = algs(idx);

end % optimal_algorithm
