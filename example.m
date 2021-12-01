function [] = example()
% Make sure the Matlab package folder `+levyarea` is
% in your current working directory or on your `path`.
% Generate a Wiener increment
m = 100;
h = 0.01;
err = 0.05;
W = sqrt(h) * randn(m,1);

% Default call chooses the best algorithm as explained
% in Section 5.1 and uses h^(3/2) as the precision.
II_a = levyarea.iterated_integrals(W,h);

% The desired precision can be provided explicitly.
II_b = levyarea.iterated_integrals(W,h,err);

% The algorithm can be chosen using a keyword argument:
% 'Fourier', 'Milstein', 'Wiktorsson' or 'MR'
II_c = levyarea.iterated_integrals(W,h,'Algorithm','Milstein');

% The error criterion can be chosen using a keyword:
% 'MaxL2' or 'FrobeniusL2'
II_d = levyarea.iterated_integrals(W,h,err,'ErrorNorm','FrobeniusL2');

% Supply square roots of eigenvalues of covariance
% operator for finite Q-Wiener processes.
q = 1./(1:m)'.^2;
QW = sqrt(h) * sqrt(q) .* randn(m,1);
II_e = levyarea.iterated_integrals(QW,h,err,'QWiener',sqrt(q));

end