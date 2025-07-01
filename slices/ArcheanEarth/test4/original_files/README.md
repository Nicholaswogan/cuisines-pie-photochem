
There are three types of reactions in `ArcheanHaze.yaml`:

## `elementary`

Suppose we have the reaction

`A + B => C + D`

The rate constant is given by `k = A*T^b*exp(-Ea/T)`, and the forward reaction rate is then given by

`Rf = [A]*[B]*k`

Where `[A]` is the number density of species A, and T is in Kelvin.

## `three-body`

Suppose we have the reaction

`A + B + M => C + M`

The rate constant is given by `k = A*T^b*exp(-Ea/T)`, and the forward reaction rate is then given by

`Rf = [A]*[B]*[M]*k`

where `[M]` is the total number density.

## `falloff`

Suppose we have the reaction

`A + B + M => C + M`

There is a low (`k0`) and high ('kinf') pressure rate constant arrhenius expression, each calculated with the usual expression (`k = A*T^b*exp(-Ea/T)`). The rate constant is then given by

`k = k0*[M]/(1 + (k0*[M]/kinf))`

and the forward reaction rate is

`Rf = [A]*[B]*k`













