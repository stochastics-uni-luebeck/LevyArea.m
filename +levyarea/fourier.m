function I = fourier(W, n)
%FOURIER Simulates Levy area using the Fourier algorithm.
%   I = FOURIER(W,N) simulates the Levy area w.r.t. the increment of 
%   a Wiener process W with variance 1.
%   Uses the Fourier algorithm with cut-off parameter N.
%
%   For a general m-dimensional Wiener increment with covariance Q=diag(q),
%   the input has to be scaled by 1/sqrt(q) and the ouput correspondingly
%   by diag(sqrt(q)) both from the left and from the right:
%       I = diag(sqrt(q)) * fourier(W./sqrt(q),n) * diag(sqrt(q)).
%
%   See also ITERATED_INTEGRALS, SIMULATE.

m = length(W);
alpha = randn(n,m);
beta = (randn(m,n) - sqrt(2).*W) ./ (1:n);
S = beta * alpha;
I = (S-S')/(2*pi);

end % fourier
