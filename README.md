# LevyArea.m
*Iterated Stochastic Integrals in Matlab*

This package implements state-of-the-art methods for the simulation of iterated stochastic integrals.
These appear e.g. in higher order algorithms for the solution of stochastic (partial) differential equations.

## Installation

Copy the folder `+levyarea` into your current working directory or into a folder on your Matlab path.

## Example

```matlab
>> h = 1/100;
>> W = sqrt(h) * randn(5,1)

W =

   -0.0875
    0.1121
    0.1414
    0.0463
    0.0907

>> II = levyarea.iteratedIntegrals(W,h, h^(3/2))

II =

   -0.0012   -0.0125   -0.0102   -0.0027   -0.0044
    0.0027    0.0013    0.0082    0.0049    0.0037
   -0.0021    0.0077    0.0050    0.0078    0.0040
   -0.0013    0.0003   -0.0012   -0.0039   -0.0029
   -0.0036    0.0065    0.0089    0.0071   -0.0009

>> 0.5*W.^2 - 0.5*h

ans =

   -0.0012
    0.0013
    0.0050
   -0.0039
   -0.0009
```
Alternatively you can import one or all functions:
```matlab
>> import levyarea.iteratedIntegrals
>> II = iteratedIntegrals(W,h, h^(3/2))

II =

   -0.0012   -0.0129   -0.0149   -0.0010   -0.0087
    0.0031    0.0013    0.0168    0.0045    0.0024
    0.0025   -0.0010    0.0050    0.0008   -0.0021
   -0.0030    0.0007    0.0058   -0.0039   -0.0036
    0.0008    0.0078    0.0149    0.0078   -0.0009
```
