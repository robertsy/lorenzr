lorenzr
=====

Personal playground for data assimilation algorithms.

Implementation of the Lorenz models:

- Lorenz63 ([Lorenz 1963](http://journals.ametsoc.org/doi/abs/10.1175/1520-0469(1963)020%3C0130:DNF%3E2.0.CO;2)), 
- Lorenz96 ([Lorenz 1998](http://journals.ametsoc.org/doi/abs/10.1175/1520-0469(1998)055%3C0399%3AOSFSWO%3E2.0.CO%3B2)), 
- Lorenz MII ([Lorenz 1998](http://journals.ametsoc.org/doi/abs/10.1175/JAS3430.1))

## Description
The core implementation is in Fortran95, with R functions provided to interface. Basic examples are given and assimilation algorithms from *assimilr* are tested.

This code was not intended as full-fledged package but rather as a quick developing platform for personal ideas concerning data assimilation algorithms. If you feel like rewriting it as a package with a better design, please contact me and I would be glad to support you if I have some time. 