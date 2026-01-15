#scratch.jl
#for quick checks and work 
# println("does this work")

include("PhysicsSimFuncs.jl")

# Solar Panel attempt
# Needed parameters, absorbed radiation, area, initial temp, surronding temp 
Tsurr = 20.0 + 273.15 # K
Ac = 1 * 2 # m^2
mdot = 0.02 # kg/s
Cp = 4186 # J/kgK (water)
Sf = 800 # W/m^2
Tif_solar = 25.0 + 273.15 # K
TiltAngle = 45 #deg 
Qf = Q_Solar(Sf, Tif_solar; Taf = Tsurr, Acf = Ac) # Returns W
println("Q to fluid: $Qf W")
Tef = Tif_solar + Qf / (mdot * Cp) # K
println("Outlet fluid temp: $(Tef - 273.15) °C")

# Tank solving 
N = 5 
Vtank = 0.3 
rho = 1000 
UA_Tank = 3 # W/K

m_node = rho * Vtank / N # kg per node

#parameters for ode 
p = (mdot = mdot, Tif = Tef, UA = UA_Tank/N, Cp = Cp, m = m_node, Ta = Tsurr)

using DifferentialEquations

T0 = fill(25.0 + 273.15, N) # Initial temp of each node in K
println("Initial Tank Energy: $(sum(T0) * m_node * Cp / 1000) kJ")
tspan = (0.0, 3600.0) # 1 hour in seconds
prob = ODEProblem(TankODE!, T0, tspan, p)
sol = solve(prob, Rodas5(), reltol=1e-6, abstol=1e-8)

TsTank = sol.u[end] # Final temperatures of each node (K) 
TotalTankEnergy = sum(TsTank) * m_node * Cp # J
println("Final Tank Energy: $(TotalTankEnergy / 1000) kJ")