#= 
PhysicsSimFuncs.jl
Jacob Child
Jan 14, 2026
Description: Functions for physics simulation
=#
using Pkg 
Pkg.activate("PassiveLogic"; shared = true)

using Infiltrator

### Solar Panel Functions ###
"""
Q_Solar(Sf, Tif; Taf = Tsurr, Acf = Ac)
Sf: Absorbed solar radiation 
Tif: Fluid inlet temperature
Taf: Ambient temperature
Acf: Area of the solar collector
"""
function Q_Solar(Sf, Tif; Taf = Tsurr, Acf = Ac)
    Tpm = Tif + 10 # Initial guess for mean plate temperature
    Ulf = UL_Solar(Tpm; Taf = Tsurr)# Overall heat loss coefficient
    Frf = FR_Solar(Ulf) # Heat removal factor
    Qu = Acf * (Sf - Ulf * (Tpm - Taf)) # Useful gain of the collector, eqn 6.9.3
    New_Tpm = Tif + (Qu / Acf) / (Frf * Ulf) * (1 - Frf) # New mean plate temperature eqn 6.9.4
    #Iterate 
    k = 0 
    while k < 100
        Tpm = New_Tpm
        Ulf = UL_Solar(Tpm; Taf = Tsurr)
        Frf = FR_Solar(Ulf)
        Qu = Acf * (Sf - Ulf * (Tpm - Taf))
        Old_Tpm = New_Tpm
        New_Tpm = Tif + (Qu / Acf) / (Frf * Ulf) * (1 - Frf)
        if abs((New_Tpm - Old_Tpm)/Old_Tpm) < 0.01
            break
        end
        k += 1
        if k == 100
            println("Warning: Max iterations reached in Q_Solar calculation")
        end
    end
    # Final calculations 
    Ulf = UL_Solar(New_Tpm; Taf = Tsurr)
    Frf = FR_Solar(Ulf)
    # Total Q to the fluid
    Qf = Acf * Frf * (Sf - Ulf*(Tif - Taf))
    # Qu_check = Acf * (Sf - Ulf * (New_Tpm - Taf))
    # if abs(Qf - Qu_check) > 1.0  # 1 W tolerance
    #     println("Warning: Energy balance mismatch: Qf=$Qf, Qu_check=$Qu_check")
    # end
     
    return Qf

end

"""
FR_Solar(Ulf; Acf = Ac, mdotf = mdot, cpf = Cp)
Ulf: Overall heat loss coefficient
Acf: Area of the solar collector
mdotf: Mass flow rate of the fluid
cpf: Specific heat capacity of the fluid
""" 
function FR_Solar(Ulf; Acf = Ac, mdotf = mdot, cpf = Cp)

    # mf = sqrt(Ulf / (kf * δf)) # eqn 6.5.4a
    # Cbf = 30 # assumed bond conductance from book and Whillier and Saluja (1965)
    # Uof = Wf * 1 / (Ulf * (Dif + (Wf-Dif)*F) + 0 + 1 / (pi * Dif * hfi) )# resistance from the absorber plate to ambient air, the '0' is the bond conductance, which I am assuming is very large, as in the book
    # Fprime = Uof / Ulf # eqn 6.5.19 in the book, ranges from 0.75 to 1.0, to simplify, I am going to arbitrarily set it to 0.875
    Fprime = 0.875 # see above for explanation, and Fig 6.5.4 in the book, this is called the collector efficiency factor
    
    cap_rate = mdotf * cpf / (Acf * Ulf * Fprime) # Dimensionless collector capacitance rate
    
    # Heat removal factor (Eq. 6.7.4)
    FRf = cap_rate * Fprime * (1 - exp(-1 / cap_rate)) # from eqn 6.7.4
    # FRf2 = (mdotf * cpf / (Acf * Ulf)) * (1 - exp(-Acf * Ulf * Fprime / (mdotf * cpf))) # eqn 6.7.4 collector capacitance rate, see section 6.7 in book
    return FRf
end

"""
UL_Solar(Tpm; Β=TiltAngle, Nf = 2, Taf = Tsurr)
Tpm: Mean plate temperature
Β: Tilt angle of the solar panel
Nf: Number of glass covers
Taf: Ambient temperature
"""
function UL_Solar(Tpm; Β=TiltAngle, Nf = 2, Taf = Tsurr)
    
    # Solve for Utop first. All below is from eqn 6.4.9
    #assuming minimal wind
    hw = 10 # W/m^2 K
    ϵg = 0.88 # emissivity of glass
    ϵp = 0.95 # emissivity of plate (assumed) 
    f = (1+0.089*hw - 0.1166*hw*ϵp)*(1+ 0.07866*Nf) 
    if Β > 70
        Β = 70
    end
    C = 520 * (1 - 0.000051*Β^2)
    ef = 0.43 * (1 - 100 / Tpm) # Tpm is the mean plate temp 
    term1 = (Nf / (C/Tpm * ((Tpm - Taf)/(Nf+f))^ef) + 1/hw)^-1
    σ = 5.67e-8 # Stefan-Boltzmann constant
    term2 = σ*(Tpm + Taf)*(Tpm^2 + Taf^2) 
    term3 = 1 / (ϵp + 0.00591*Nf*hw)
    term4 = (2*Nf + f - 1 + 0.133 * ϵp) / ϵg
    Utop = term1 + term2 / (term3 + term4 - Nf)
    Ubot = 0.9 #arbitrary, from book example should be calculated as k/L for back wall insulation material

    return Utop + Ubot
end

### Tank Functions ### 
"""
function Fc(i, Tfi, Ts)
i: Index of the current tank node
Tfi: Inlet fluid temperature
Ts: Array of tank temperatures
"""
@inline function Fc(i, Tfi, Ts)
    N = length(Ts)

    if i == 1 && Tfi > Ts[1]
        return 1.0
    elseif i >= 2 && Ts[i-1] >= Tfi > Ts[i]
        return 1.0
    elseif i == 0 || i == N + 1
        return 0.0
    else 
        return 0.0
    end
end

# For now I am going to assume no load, so I will ignore the load terms 
"""
function mixed_flow_mdots(Fcs, mdot) #eqn 8.4.3
Fcs: Array of flow control factors
mdot: Mass flow rate of the fluid
"""
function mixed_flow_mdots(Fcs, mdot) #eqn 8.4.3
    N = length(Fcs)
    mdot_m = zeros(N+1)
    for i in 2:N # indexing seems different, but [1] = 0 and [N+1] = 0 as they should 
        mdot_m[i] = mdot * sum(Fcs[1:i-1])
    end
    return mdot_m
end 


"""
ode form:
function f!(dT, T, p, t)
    where dT is the derivative of T with respect to time and p are parameters
    the function updates dT in place ie dT[i] = ...
end
"""
function TankODE!(dT, Ts, p, t)
    mdot, Tif, UA, Cp, m, Ta = p 
    N = length(Ts)

    Fcs = [Fc(i, Tif, Ts) for i in 1:N] # Flow control function 
    

    mdot_ms = mixed_flow_mdots(Fcs, mdot) # Mixed flow mdots

    for i in 1:N 
        Losses = (UA / Cp) * (Ta - Ts[i])
        InletTerm = Fcs[i] * mdot * (Tif - Ts[i])
        #mixing 
        MixingTerm = 0.0
        if mdot_ms[i] > 0 && i > 1
            MixingTerm = mdot_ms[i] * (Ts[i-1] - Ts[i])
        elseif mdot_ms[i+1] < 0 && i < N
            MixingTerm = mdot_ms[i+1] * (Ts[i] - Ts[i+1])
        end
        dT[i] = (Losses + InletTerm + MixingTerm) / m # this assumes every section is the same mass
    end
end