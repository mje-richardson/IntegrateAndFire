# IntegrateAndFire
Julia code for rate and densities for the leaky and exponential integrate-and-fire neurons

Version date: 12/01/2024 using Julia 1.9   
Provides rates and probability densities for the LIF and EIF models and includes steady state and current (mu) modulation.  
For the LIF: tau×dvdt=mu - v + sig×sqrt(2×tau)×xi(t)  
For the EIF: tau×dvdt=mu - v + dT×exp((v-vT)/dT) + sig×sqrt(2×tau)×xi(t)  
where xi(t) is gaussian white noise with unit amplitude Dirac-delta autocov.  

Note: the code uses a second-order exponential scheme that generalises the algorithms in  
Richardson MJE. Phys Rev E 76 021919 (2007)
