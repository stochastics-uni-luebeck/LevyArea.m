function G = mronroe(W, n)
%MRONROE Simulates Levy area using the Mrongowius-Roessler algorithm.
%   I = MRONROE(W,N) simulates the Levy area w.r.t. the increment of 
%   a Wiener process W with variance 1.
%   Uses the Mrongowius-Roessler algorithm with cut-off parameter N.
%
%   For a general m-dimensional Wiener increment with covariance Q=diag(q),
%   the input has to be scaled by 1/sqrt(q) and the ouput correspondingly
%   by diag(sqrt(q)) both from the left and from the right:
%       I = diag(sqrt(q)) * mronroe(W./sqrt(q),n) * diag(sqrt(q)).
%
%   See also ITERATED_INTEGRALS, SIMULATE.

m = length(W);
alpha = randn(n,m);
beta = (randn(m,n) - sqrt(2).*W) ./ (1:n);
S = beta * alpha;
M = randn(1,m);
G = zeros(m);
G(tril(true(m),-1)) = randn(m*(m-1)/2,1);
S = S + sqrt(2*psi(1,n+1)) * (W.*M + G);
G = (S-S')/(2*pi);

end % mronroe
