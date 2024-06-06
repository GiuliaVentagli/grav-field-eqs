# Manipulating and integrating the field equations for GR or alternative theories of gravity with Julia Symbolics.

This repository contains the code for writing in component notations the field equations for either General Relativity or a generic scalar-tensor theory starting from tensorial quantities in Julia with the use of the Symbolics package. The equations can also be manipulated and put in the form $x'=F(x,r)$, and converted automatically into Julia code to define a solvable ODEs system. The codes assumes $c=1$

## Packages used

The current version of Julia used is 1.10.3. The versions of packages used in this repository are as follow:
* Symbolics v5.28.0
* LinearAlgebra 
* DifferentialEquations v7.13.0
* Plots v1.40.4

## GravityTensors.jl

All quantities, apart from the needed variables, the metric and the inverse metric, are defined in the module "GravityTensors.jl". Here is a list of the functions currently contained in the module:
* derMet(*Met*): computes the partial derivative of a given metric *Met*, $\partial_\mu g_{\mu\nu}$.
* christoffel(*Met*,*InvMet*): computes the Christoffel symbol given a metric *Met* and its inverse *InvMet*, $\Gamma^\kappa_{\mu\nu}$.
* derChris(*chris*): computes the partial derivative of the Christoffel symbol *chris*, $\partial_\nu \Gamma^\kappa_{\mu\rho}$.
* riemann(*chris*): computes the Riemann tensor given the Christoffel symbol *chris*, $R^\rho_{\mu\nu\sigma}$.
* ricciTens(*riem*): computes the Ricci tensor given the Riemann tensor *riem*, $R_{\mu\nu}$.
* ricciScal(*Met*,*ricciT*): computes the Ricci scalar given the metric *Met* and the Ricci tensor *ricciT*, $R=g^{\mu\nu}R_{\mu\nu}$.
* gaussBonnet(*Met*,*InvMet*,*riem*,*ricciT*,*scalar*): computes the Gauss-Bonnet invariant given the metric *Met*, the inverse metric *InvMet*, the Riemann tensor *riem*, the Ricci tensor *ricciT* and the Ricci scalar *scalar*, $\mathcal{G}=R^2-4R^{\mu\nu}R_{\mu\nu}+R^{\mu\nu\rho\sigma}R_{\mu\nu\rho\sigma}$.
* stressEnergy(*Met*,*vel*, $\epsilon$ ,*p*): computes the perfect fluid stress-energy tensor given the metric *Met*, the 4-velocity *vel*, the energy density $\epsilon$ and pressure *p*, $T^{\mu\nu}=(\epsilon + p)u^\mu u^\nu+p g^{\mu\nu}$.
* stressEnergydd(*Met*,*T*): computes the covariant perfect fluid stress-energy tensor given the metric *Met* and the controvariant stress-energy tensor *T*, $T_{\mu\nu}$.
* DerStressEn(*chris*,*T*): computes the conservation of the stress-energy tensor given the Christoffel symbol *chris* and the controvariant stress-energy tensor *T*, $\nabla_\mu T^{\mu\nu}$.
* derDerScal(*chris*, $\phi$): computes the double covariant derivative of a scalar field given the Christoffel symbol *chris* and the scalar field $\phi$, $\nabla_\mu\nabla_\nu\phi$.
* boxScal(*InvMet*,*chris*, $\phi$): computes the d'Alembertian of a scalar field given the inverse of the metric *InvMet*, the Christoffel symbol *chris* and the scalar field $\phi$, $\Box\phi$.
* dScaldScal($\phi$): computes the product of two derivatives of a scalar field given the scalar field $\phi$, $\nabla_\mu\phi\nabla_\nu\phi=\partial_\mu\phi\partial_\nu\phi$.
* dScaldScalContr(*InvMet*, $\phi$): computes the contracted product of two derivatives of a scalar field given the inverse of the metric *InvMet* and the scalar field $\phi$, $\nabla_\mu\phi\nabla^\mu\phi=\partial_\mu\phi\partial^\mu\phi$.

## NSexample.ipynb

This is a practical example that uses "GravityTensors.jl". The code solves the field equations for a static spherically symmetric neutron star, described by an equation of state given by a polytrope. After computing the needed quantities, we define the field equations, we manipulate them and we create numerically-usable functions in Julia language. We then solved the system and produce the pressure vs radial coordinate plot, as an example.
