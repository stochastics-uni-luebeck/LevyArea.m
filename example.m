function [] = example()
% Make sure the Matlab files are in your current working
% directory or in your `path`.
% Generate a Wiener increment
m = 500;
h = 0.001;
eps = 0.00001; % = h^(3/2)
W = sqrt(h) * randn(m,1);

% Default call chooses the best algorithm as explained
% in Section 5.1 and uses h^(3/2) as the precision.
II = levyarea.iteratedIntegrals(W,h);
% The desired precision can be provided explicitly.
II = levyarea.iteratedIntegrals(W,h,eps);
% The algorithm can be chosen using a keyword argument:
% 'Fourier', 'Milstein', 'Wiktorsson' or 'MR'
II = levyarea.iteratedIntegrals(W,h,eps, 'Algorithm', 'Fourier');
% Supply eigenvalues of square root of covariance matrix
% for finite Q-Wiener processes.
q = 1./(1:m)'.^2;
QW = sqrt(h) * sqrt(q) .* randn(m,1);
II = levyarea.iteratedIntegrals(QW,h,eps,'QWiener',sqrt(q));

end