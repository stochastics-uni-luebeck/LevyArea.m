function I = iteratedIntegrals(W, h, varargin)
%ITERATEDINTEGRALS Simulates twice iterated stochastic integrals.
%   I = ITERATEDINTEGRALS(W,H,EPS) simulates the twice iterated stochastic
%   integral w.r.t. the increment of a Wiener process W over stepsize H.
%   The algorithm guarantees that the error is less than or equal to EPS in
%   the specified error norm.
%   For SDEs the Maximum-L2 error estimate is used, whereas for SPDEs 
%   the Frobenius-L2 error estimate is used.
%
%   If the integrals are to be used in a SDE solver scheme then EPS should
%   typically be chosen as H^(p+0.5) where p is the order of the SDE scheme.
%   
%   I = ITERATEDINTEGRALS(W,H,EPS,PARAM1,VAL1,PARAM2,VAL2,...) performs
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
%

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
addParameter(ip,'MemUse',1,@(x) isscalar(x) && isnumeric(x) && ...
    0<=x && x<=1);
addParameter(ip,'QWiener',ones(m,1),@(x) isnumeric(x) && ...
    length(x)==m && iscolumn(x));
addParameter(ip,'ErrorNorm',"auto",@(x) isstring(x) || ischar(x))
parse(ip,W,h,varargin{:});

eps = ip.Results.Error;
if ip.Results.Algorithm == "auto"
    alg = levyarea.optimalAlgorithm(m,h,eps);
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

I = levyarea.simulate(W,h,eps,ip.Results.ItoCorrection,alg,...
    ip.Results.MemUse,ip.Results.QWiener,err_norm);

end % iteratedIntegrals