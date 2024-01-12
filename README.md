# IntegrateAndFire
Julia code for rate and densities for the leaky and exponential integrate-and-fire neurons

Version date: 12/01/2024 using Julia 1.9
Provides rates and probability densities for the LIF and EIF models and includes steady state and current (mu) modulation.
For the LIF: tau*dvdt=mu-v+2*sig*sqrt(tau)*xi(t)
For the EIF: tau*dvdt=mu-v+dT*exp((v-vT)/dT)+2*sig*sqrt(tau)*xi(t)
where xi(t) is gaussian white noise with unit amplitude Dirac-delta autocov.

Note: the code uses a 2nd-order exponential scheme tgat generalises the algorithms in Richardson MJE. Phys Rev E 76 021919 (2007)
