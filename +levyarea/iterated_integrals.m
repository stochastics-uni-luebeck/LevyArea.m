function I = iterated_integrals(W, h, varargin)
%ITERATED_INTEGRALS Simulates twice iterated stochastic integrals.
%   I = ITERATED_INTEGRALS(W,H) simulates the twice iterated stochastic
%   integral w.r.t. the increment of a Wiener process W over stepsize H.
%   The algorithm guarantees that the error is less than or equal to
%   H^(3/2) in the Max-L2 error norm.
%
%   I = ITERATED_INTEGRALS(W,H,Err) simulates the twice iterated stochastic
%   integral w.r.t. the increment of a Wiener process W over stepsize H.
%   The algorithm guarantees that the error is less than or equal to Err in
%   the Max-L2 error norm.
%
%   If the integrals are to be used in a SDE solver scheme then ERR should
%   typically be chosen as H^(p+0.5) where p is the order of the SDE scheme.
%   
%   I = ITERATED_INTEGRALS(___,Name,Value) performs
%   the simulation with the specified values of optional parameters.
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
%       Possible values are 'Auto', 'Fourier', 'Milstein', 'Wiktorsson'
%       and 'MronRoe'. 'Auto' chooses the algorithm according to
%       `levyarea.optimal_algorithm`. 
%       The default value is 'Auto'.
%
%   'QWiener', square roots of the eigenvalues of the covariance operator
%              of a Q-Wiener process with Q=diag(QWiener)^2
%   
%       When 'QWiener' is supplied, W is assumed to be a finite-dimensional
%       approximation of an infinite Q-Wiener process with covariance
%       matrix Q=diag(QWiener)^2. 
%       The default value is a vector of ones.
%
%   'ErrorNorm', which norm to use for the error estimate
%
%       Possible values are 'Auto', 'MaxL2' and 'FrobeniusL2'. For SDEs usually
%       the Maximum-L2 error estimate is needed, whereas for SPDEs the 
%       Frobenius-L2 error estimate is needed. 'Auto' chooses 'MaxL2' if
%       'q_12' is a vector of ones and 'FrobeniusL2' otherwise.
%       The default value is 'Auto'.
%       
%   See also OPTIMAL_ALGORITHM.

% dimension of the Wiener process
m = length(W);

% check input arguments
narginchk(2,inf);
ip = inputParser;
addRequired(ip,'WienerIncrement',@(x) isnumeric(x) && iscolumn(x));
addRequired(ip,'StepSize',@(x) (x>0) && isnumeric(x) && isscalar(x));
addOptional(ip,'Error',h^(3/2),@(x) (x>0) && isnumeric(x) && isscalar(x));
addParameter(ip,'ItoCorrection',true,@(x) isscalar(x) && ...
    (islogical(x) || (isnumeric(x) && (x == 0 || x == 1))));
addParameter(ip,'Algorithm',"auto",...
    @(x) isstring(x) || ischar(x));
addParameter(ip,'QWiener',ones(m,1),@(x) isnumeric(x) && ...
    length(x)==m && iscolumn(x));
addParameter(ip,'ErrorNorm',"auto",@(x) isstring(x) || ischar(x))
parse(ip,W,h,varargin{:});

err = ip.Results.Error;
if ip.Results.Algorithm == "auto"
    alg = levyarea.optimal_algorithm(m,h,err);
else
    alg = ip.Results.Algorithm;
end
if ip.Results.ErrorNorm == "auto"
    if all(ip.Results.QWiener == 1)
        err_norm = "maxl2";
    else
        err_norm = "frobeniusl2";
    end
else
    err_norm = ip.Results.ErrorNorm;
end

I = levyarea.simulate(W,h,err,ip.Results.ItoCorrection,alg,...
    ip.Results.QWiener,err_norm);

end % iterated_integrals