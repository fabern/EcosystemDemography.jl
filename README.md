# EcosystemDemography.jl

A barebone implementation of a PDE representing a simple Forest Ecosystem Demography Model. Terminology follows the definitions in Lehmann et al. (2015, http://dx.doi.org/10.1016/j.ecolmodel.2015.01.013).

To learn more about ecosystem demography models in general refer to the description in the book 'Climate Change and Terrestrial Ecosystem Modeling' by Bonan (2019, https://doi.org/10.1017/9781107339217).


`01_intro_MOL.jl` is a short intor how to solve PDEs in Julia with the package `MethodsOfLines.jl`.

`02_intro_ED.jl` solves the model from [Lehmann et al. (2015)](http://dx.doi.org/10.1016/j.ecolmodel.2015.01.013), i.e.


## Ecosystem dynamics model from [Lehmann et al. (2015)](http://dx.doi.org/10.1016/j.ecolmodel.2015.01.013), 
### Base partial differential equation
Ecosystem dynamics model from [Lehmann et al. (2015)](http://dx.doi.org/10.1016/j.ecolmodel.2015.01.013), i.e.
```math
\frac{\partial n_i(t,x)}{\partial t} = -\frac{\partial g_i(t,x)n_i(t,x)}{\partial x} - m_i(t,x)n_i(t,x),
```

where $n_i(t,x)$ is the number of trees per area of species $i$ (expressed as density in dbh-space, i.e with units of cm⁻¹ m⁻²) as a function of time $t$ (years) and diameter at breast height (DBH) $x$ (cm).

The time evolution of the species size distribution is then determined by the species-specific growth and mortality rate functions $g_i(t,x)$ (cm year⁻¹) and $m_i(t,x)$ (year⁻¹) as well as by the recruitment function $r_i(t)$ (m⁻² year⁻¹).

Boundary conditions introduce the recruitment of small trees of DBH size $x_0$, i.e. 
```math
n_i(t,x_0)g_i(t,x_0) = r_i(t).
```

Initial conditions are specified as $n_i(0,x) = 0$ (cm⁻¹ m⁻²).

### Inventory parameters
Inventory parameters can be derived as:
- number of trees (m⁻²) in diameter class [ $x_a$, $x_b$ ] (as the integral of the density): 
  
  $N_i(t, x_a, x_b) = \int_{x_a}^{x_b} n_i(t, x) dx,$
- total number of trees above $x_0$: 
  
  $N_i(x_0, \infty)$ (m⁻²)
- cumulative basal area (modified from Lehmann 2015) of trees larger than $x$:
  
  $BA_i(t, x) = \int_{x}^{\infty} \frac{\pi}{4} y^2 n_i(t,y) dy $
- total basal area of trees above $x_0$:

  $BA_i(t, x_0)$
- and of all species 
  
  $BA(t,x) = \sum_i BA_i(t,x)$
- total aboveground biomass
  
  $BM = TODO$


### Models for growth, mortality, and recruitment
#### Constant in time and constant with $x$
With constant growth $g_i$, mortality $m_i$, and recruitment $r_i$ the steady-state distribution 
becomes exponential, i.e.
```math
\frac{n_i(x)}{N_i} = \frac{m_i}{g_i} \exp\left(- \frac{m_i}{g_i}x\right),
```
where $N_i$ is determined by the recruitment and mortality rates $N_i = \frac{r_i}{m_i}$.

TODO: we can use this to test the implementation.

#### Constant in time and varying with $x$
According to Lehmann et al. (2015) the steady state solution with of time-constant $g_i(x)$, $m_i(x)$, and $r_i(x_0)$ is:
```math
n_i(x) = \frac{r_i}{g_i(x_0)} \exp\left(-\int_{x_0}^{x}\frac{m_i(y)}{g_i(y)}dy\right).
```

TODO: we can use this to test the implementation.

#### Lehmann et al. (2015)
In Lehmann et al. (2015) they model mortality and recruitment to be constant, but make the growth dependent on the competition for light by trees in larger size classes.

This can be done by specifying:
```math
g_i(t,x) = g_{i0} \cdot max(0, (1-\frac{BA(t,x)}{\beta_i})),
```

where $BA(t,x)$ refers to the cumulative basal area of all species above the DBH $x$.

Then they solve this numerically with an upwind scheme.

TODO: e.g when using Julia we can do this with MethodOfLines.jl capable of solving [PIDEs](https://sciml.github.io/MethodOfLines.jl/dev/tutorials/PIDE/#integral).

NOTES: there exist tailor-made space discretizations for the general type of these equations (McKendrick-von Foerster Partial Differential Equation (PDE)). 
E.g. see the methods implemented by [libpspm](https://jaideep777.github.io/libpspm/articles/what_is_pspm.html). 
(NOTE however, that this appears to support PDEs but not PIDEs, i.e. interaction between size classes or between species is not implemented. 
NOTE: this might be incorrect. There appears to be support for integrals from `x` to `xmax` by using [`integrate_wudx_above()`](https://jaideep777.github.io/libpspm/articles/size_integral.html#computing-partial-state-integrals).)

#### BiomeE
A simplified version of [BiomeE](https://github.com/geco-bern/rsofun) could be implemented by making a PDE in dimensions of DBH $x$ and height $y$, and then reducing the growth by the cumulative LAI of all trees above height $y$. TODO....

#### P-Model dependent growth
TODO....


