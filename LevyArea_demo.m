% Example script for LevyArea.m

% Make sure the MATLAB toolbox folder `LevyArea.m` is on your `path`
LevyArea_setup

% Generate a Wiener increment
m = 100;
h = 0.01;
err = 0.05;
W = sqrt(h) * randn(m,1);

% Default call chooses the best algorithm as explained
% in Section 5.1 and uses h^(3/2) as the precision.
II = levyarea.iterated_integrals(W,h);

% The desired precision can be provided explicitly.
II = levyarea.iterated_integrals(W,h,err);

% The algorithm can be chosen using a keyword argument:
% 'Fourier', 'Milstein', 'Wiktorsson' or 'MronRoe'
II = levyarea.iterated_integrals(W,h,'Algorithm','Milstein');

% The error criterion can be chosen using a keyword:
% 'MaxL2' or 'FrobeniusL2'
II = levyarea.iterated_integrals(W,h,err,'ErrorNorm','FrobeniusL2');

% Supply square roots of eigenvalues of covariance
% operator for finite Q-Wiener processes.
q = 1./(1:m)'.^2;
QW = sqrt(h) * sqrt(q) .* randn(m,1);
IIQ = levyarea.iterated_integrals(QW,h,err,'QWiener',sqrt(q));
