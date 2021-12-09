function I = simulate(W, h, err, ito, alg, q_12, err_norm)
%SIMULATE The core routine to simulate iterated integrals.
%   I = SIMULATE(W,H,Err,Ito,Alg,QHalf,ErrNorm) calculates the correct
%   cut-off parameter for the chosen algorithm using the given information
%   and then calls that algorithm. This function also handles the scaling
%   of the iterated integrals, i.e. W can be a Q-Wiener increment over
%   stepsize H. Return value is the matrix of iterated integrals (not the
%   Levy areas).
%
%   This function does no special parameter handling, i.e. all arguments
%   are required and positional arguments.
%
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
        I = diag(sqcov)*levyarea.fourier(W./sqcov,n)*diag(sqcov);
    case "milstein"
        n = ceil( 0.5*(norm_coeff*h/(pi*err))^2 );
        I = diag(sqcov)*levyarea.milstein(W./sqcov,n)*diag(sqcov);
    case "wiktorsson"
        n = ceil( sqrt(5*m/12) * norm_coeff*h/(pi*err) );
        I = diag(sqcov)*levyarea.wiktorsson(W./sqcov,n)*diag(sqcov);
    case "mronroe"
        n = ceil( sqrt(m/12) * norm_coeff*h/(pi*err) );
        I = diag(sqcov)*levyarea.mronroe(W./sqcov,n)*diag(sqcov);
    otherwise
        error("Unknown algorithm for Q-Wiener processes: " + ...
            ip.Results.Algorithm + ...
            ". Possible choices are: Fourier, Milstein, Wiktorsson, MronRoe.");
end

% calculate iterated integrals
I = I + W.*W'./2;

% optionally add Ito correction
if ito
    I = I - diag(q_12).^2.*h./2;
end

end % simulate
