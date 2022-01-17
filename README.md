# LevyArea.m
*Iterated Stochastic Integrals in Matlab*

This package implements state-of-the-art methods for the simulation of iterated stochastic integrals.
These appear e.g. in higher order algorithms for the solution of stochastic (partial) differential equations.

## Installation

Install the Matlab toolbox file (`LevyArea.mltbx`) from the [Releases](https://github.com/stochastics-uni-luebeck/LevyArea.m/releases) page or through the Add-On Explorer.
Alternatively you can copy the folder `+levyarea` into your current working directory or into a folder on your Matlab path.

## Example

The main function of the toolbox is the function
`iterated_integrals`. It can be called by prepending the
package name `levyarea.iterated_integrals`. However, since
this may be cumbersome one can import the used functions once by
```matlab
>> import levyarea.iterated_integrals
>> import levyarea.optimal_algorithm
```
and then one can omit the package name by simply calling
`iterated_integrals` or `optimal_algorithm`,
respectively.
In the following, we assume that the two functions are imported, so that we can always call them directly without the package name.


First we generate a Wiener increment:
```matlab
>> m = 5;	% dimension of Wiener process
>> h = 0.01;	% step size or length of time interval
>> err = 0.05;	% error bound
>> W = sqrt(h) * randn(m,1);	% increment of Wiener process
```
Here, $W$ is the $m$-dimensional vector of increments of the driving
Wiener process on some time interval of length $h$.

The default call uses h^(3/2) as the precision and chooses the best algorithm automatically:
```matlab
>> II = iterated_integrals(W,h)
```
If not stated otherwise, the default error criterion is the $\max,L^2$-error
and the function returns the $m \times m$ matrix `II` containing a realisation
of the approximate iterated stochastic integrals that correspond
to the given increment $W$.

The desired precision can be optionally provided
using a third positional argument:
```matlab
>> II = iterated_integrals(W,h,err)
```
Again, the software package automatically chooses the optimal
algorithm.

To determine which algorithm is chosen by the package without simulating any iterated
stochastic integrals yet, the function `optimal_algorithm` can
be used. The arguments to this function are the dimension of the Wiener
process, the step size and the desired precision:
```matlab
>> alg = optimal_algorithm(m,h,err); % output: 'Fourier'
```

It is also possible to choose the algorithm directly using a key-value pair.
The value can be
one of `'Fourier'`, `'Milstein'`, `'Wiktorsson'`
and `'MronRoe'`:
```matlab
>> II = iterated_integrals(W,h,'Algorithm','Milstein')
```

The desired norm for the prescribed error bound can also be
selected using a key-value pair. The accepted values are
`'MaxL2'` and `'FrobeniusL2'` for the $\max,L^2$- and $\mathrm{F},L^2$-norm, respectively:
```matlab
>> II = iterated_integrals(W,h,err,'ErrorNorm','FrobeniusL2')
```

The simulation of numerical solutions to SPDEs often requires
iterated stochastic integrals based on $Q$-Wiener processes.
In that case, the square roots of the eigenvalues of the associated
covariance operator need to be provided:
```matlab
>> q = 1./(1:m)'.^2; % eigenvalues of covariance operator
>> QW = sqrt(h) * sqrt(q) .* randn(m,1); % Q-Wiener increment
>> IIQ = iterated_integrals(QW,h,err,'QWiener',sqrt(q))
```
In this case, the function utilizes a
scaling of the iterated stochastic integrals and also adjusts the error
estimates appropriately such that the error bound holds w.r.t.\ the
iterated stochastic integrals $\mathcal{I}^{Q}(h)$ based on the
$Q$-Wiener process. Here the error norm defaults to the $\mathrm{F},L^2$-error.

Note that all discussed keyword arguments (key-value pairs) are
optional and can be combined as needed. Additional information
can be found using the `help` function:

```matlab
>> help iterated_integrals
>> help levyarea.optimal_algorithm
```

## Related Packages

A Julia version of this package is also available under [LevyArea.jl](https://github.com/stochastics-uni-luebeck/LevyArea.jl).