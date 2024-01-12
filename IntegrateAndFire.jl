################################################################################
################################################################################
################################################################################
################################################################################
# IntegrateAndFire.jl 
# Version date: 12/01/2024 using Julia 1.9
# Published under GNU General Public License v3.0
#
# Provides rates and probability densities for the LIF and EIF models and
# includes steady state and current (mu) modulation.
# For the LIF: tau*dvdt=mu-v+2*sig*sqrt(tau)*xi(t)
# For the EIF: tau*dvdt=mu-v+dT*exp((v-vT)/dT)+2*sig*sqrt(tau)*xi(t)
# where xi(t) is gaussian white noise with unit amplitude Dirac-delta autocov.
#
# Note: the code uses a 2nd-order exponential scheme tgat generalises 
# the algorithms in Richardson MJE. Phys Rev E 76 021919 (2007)
################################################################################
################################################################################
################################################################################
################################################################################

module IntegrateAndFire
export SteadyStateLIF, RateModLIF
export SteadyStateEIF, RateModEIF

################################################################################
################################################################################
################################################################################
################################################################################
# LIF code
################################################################################
################################################################################
################################################################################
################################################################################
 
################################################################################
################################################################################
# LIF steady state
################################################################################
################################################################################
function SteadyStateLIF(tau,mu,sig,vlo,vre,vth,dv)

    # set up voltage range and find reset points
    v=[collect(vlo:dv:vre);collect(vre:dv:vth)];
    n=length(v); krem,krep=findall(v.==vre)
    
    # create useful arrays 
    phi=(v.-mu).^2/(2*sig^2); dphi=(v.-mu)/sig^2; 
    a=exp.([0;diff(phi)]); 
    I=(exp.(dphi*dv) .-1)./dphi; I[findall(dphi.==0)].=dv
    
    # arrays with boundary conditions at vth
    p0,j0=zeros(n),ones(n)

    for k=n:-1:2
        
        # quantities needed for integration
        hk=tau*j0[k]/sig^2
        bk=hk*I[k-1]
        
        # the integration
        p0[k-1]=a[k]*p0[k] + bk*(k!==krep)
        j0[k-1]=j0[k] - (k==krep)
        
    end
    
    dp0dv=-dphi.*p0-tau*j0/sig^2
    r0=1/(sum(p0[1:n]*dv) - p0[krep]*dv);
    
    return r0,v,r0*p0,r0*j0,r0*dp0dv
end

################################################################################
################################################################################
# LIF rate modulation
################################################################################
################################################################################
function RateModLIF(tau,mu,sig,vlo,vre,vth,w,dv)
    
    # get the steady state
    r0,v,P0,J0,dP0dv=SteadyStateLIF(tau,mu,sig,vlo,vre,vth,dv)
    n=length(v); krem,krep=findall(v.==vre)
    
    # useful arrays
    phi=(v.-mu).^2/(2*sig^2); dphi=(v.-mu)/sig^2; 
    a=exp.([0;diff(phi)]); 
    I=(exp.(dphi*dv) .-1)./dphi; L=(I .-dv)./dphi;
    I[findall(dphi.==0)].=dv; L[findall(dphi.==0)].=dv^2/2
    
    # inhomogeneus solution for the modulated rate
    jr=ones(Complex,n); pr=zeros(Complex,n);

    # inhomogeneous solution for the driving modulated mu
    ju=zeros(Complex,n); pu=zeros(Complex,n);

    for k=n:-1:2
        
        ################################################################
        # Solution for the rate part
        ################################################################
        
        # gradients
        dprdvk=-dphi[k]*pr[k]-(tau/sig^2)*jr[k]
        djrdvk=-1im*w*pr[k]
        # prefactors
        hk=tau*jr[k]/sig^2; 
        dhk=tau*djrdvk/sig^2
        bk=I[k-1]*hk - L[k-1]*dhk
        ck=1im*w*(dv*pr[k]-(dv^2/2)*dprdvk)
        # iteration
        pr[k-1]=a[k]*pr[k] + bk*(k!==krep)
        jr[k-1]=jr[k] + ck*(k!==krep) - (k==krep)
        
        ################################################################
        # Solution for the mu part
        ################################################################

        # gradients needed
        dpudvk=-dphi[k]*pu[k] - (tau/sig^2)*ju[k] + P0[k]/sig^2
        djudvk=-1im*w*pu[k]
        # prefactors        
        hk=tau*ju[k]/sig^2 - P0[k]/sig^2; 
        dhk=tau*djudvk/sig^2 - dP0dv[k]/sig^2
        bk=I[k-1]*hk - L[k-1]*dhk
        ck=1im*w*(dv*pu[k]-(dv^2/2)*dpudvk)
        # iteration        
        pu[k-1]=a[k]*pu[k] + bk*(k!==krep)
        ju[k-1]=ju[k] + ck*(k!==krep)
    
    end
       
    # find the modulated rate
    rh=-ju[1]/jr[1]; # because Jh[vlo]=0 => ju[1] + r1*jr[1]=0;
    ph=pr*rh .+ pu; 
    jh=jr*rh .+ ju;

    return r0,rh,v,ph,jh

end

################################################################################
################################################################################
################################################################################
################################################################################
# EIF code
################################################################################
################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
# EIF steady state
################################################################################
################################################################################
function SteadyStateEIF(tau,mu,dT,vT,sig,vlo,vre,vth,dv)
    
    # set up voltage range and find reset points
    v=[collect(vlo:dv:vre);collect(vre:dv:vth)]; 
    n=length(v); krem,krep=findall(v.==vre)
    
    # create useful arrays  
    phi=((v .-mu).^2/2 - dT^2*exp.((v .-vT)/dT))/sig^2
    dphi=(v .- mu .- dT*exp.((v .-vT)/dT))/sig^2 ;
    a=exp.([0;diff(phi)]); 
    I=(exp.(dphi*dv) .-1)./dphi; I[findall(dphi.==0)].=dv
    
    # arrays with boundary conditions at vth
    p0,j0=zeros(n),ones(n); 
    p0[n]=-tau/(sig^2*dphi[n])  # an initial estimate
    
    for k=n:-1:2
        
        # quantities needed for integration
        hk=tau*j0[k]/sig^2
        bk=hk*I[k-1]
        
        # the integration
        p0[k-1]=a[k]*p0[k] + bk*(k!==krep)
        j0[k-1]=j0[k] - (k==krep)
        
    end

    dp0dv=-dphi.*p0-tau*j0/sig^2
    r0=1/(sum(p0[1:n]*dv) - p0[krep]*dv);
    
    return r0,v,r0*p0,r0*j0,r0*dp0dv
end

################################################################################
################################################################################
# EIF rate modulation
################################################################################
################################################################################
function RateModEIF(tau,mu,dT,vT,sig,vlo,vre,vth,w,dv)
    
    # get the steady state
    r0,v,P0,J0,dP0dv=SteadyStateEIF(tau,mu,dT,vT,sig,vlo,vre,vth,dv)
    n=length(v); krem,krep=findall(v.==vre)
    
    # useful arrays
    phi=((v .-mu).^2/2 - dT^2*exp.((v .-vT)/dT))/sig^2
    dphi=(v .- mu .- dT*exp.((v .-vT)/dT))/sig^2 ;
    a=exp.([0;diff(phi)]); 
    I=(exp.(dphi*dv) .-1)./dphi; L=(I .-dv)./dphi;
    I[findall(dphi.==0)].=dv; L[findall(dphi.==0)].=dv^2/2
    
    # inhomogeneus solution for the modulated rate
    jr=ones(Complex,n); pr=zeros(Complex,n);
    pr[n]=-tau/(sig^2*dphi[n])  # an initial estimate

    # inhomogeneous solution for the driving modulated mu
    ju=zeros(Complex,n); pu=zeros(Complex,n);

    for k=n:-1:2
    
        ################################################################
        # Solution for the rate part
        ################################################################
        
        # gradients 
        dprdvk=-dphi[k]*pr[k]-(tau/sig^2)*jr[k]
        djrdvk=-1im*w*pr[k]
        # prefactors
        hk=tau*jr[k]/sig^2
        dhk=tau*djrdvk/sig^2
        bk=I[k-1]*hk - L[k-1]*dhk
        ck=1im*w*(dv*pr[k]-(dv^2/2)*dprdvk)
        # iteration
        pr[k-1]=a[k]*pr[k] + bk*(k!==krep)
        jr[k-1]=jr[k] + ck*(k!==krep) - (k==krep)
        
        ################################################################
        # solution for the mu mod part
        ################################################################
        
        # gradients needed
        dpudvk=-dphi[k]*pu[k] - (tau/sig^2)*ju[k] + P0[k]/sig^2
        djudvk=-1im*w*pu[k]
        # prefactors
        hk=tau*ju[k]/sig^2-P0[k]/sig^2
        dhk=tau*djudvk/sig^2-dP0dv[k]/sig^2
        bk=I[k-1]*hk - L[k-1]*dhk
        ck=1im*w*(dv*pu[k]-(dv^2/2)*dpudvk)
        # iteration
        pu[k-1]=a[k]*pu[k] + bk*(k!==krep)
        ju[k-1]=ju[k] + ck*(k!==krep)
    
    end
        
    # find the modulated rate
    rh=-ju[1]/jr[1]; # because Jh[vlo]=0 => ju[1] + r1*jr[1]=0;
    ph=pr*rh .+ pu; 
    jh=jr*rh .+ ju;

    return r0,rh,v,ph,jh
    
end

end