I am to write a simple software simulation of the following system. 
Minimum Requirements
1. The system should simulate the heat transfer from a solar panel to a storage tank
2. Any coding language
3. Evaluated on thermodynamic correctness, code approach, and results

System: Environment -> Solar Panel/Heat exchanger pipe system -> pump -> top of storage tank -> cooler water from bottom of storage tank -> solar panel :||

To start I might follow l[this](https://allendowney.github.io/ModSimPy/chap01.html) course.

# Equations
[Solar Engineering of Thermal Processes](https://lib.byu.edu/search/byu/record/cram.126.ocn836402985.1) chp 12, 
# CV1 (Solar Panel)
$Q_u = A_c F_R [S - U_L ( T_i - T_a)]$ From [Penn State](https://courses.ems.psu.edu/eme811/node/730), or eqn 6.7.6 from the book above
Where $F_R$ is the heat removal factor found in section 6.7 of the book, needs collector efficiency factor $F'$
$F'$ can be found from 6.5.19 $F' = U_o / U_L$ where $U_o$ is the resistance from the absorber plate to ambient air (look up how to find this)
$A_c$ Area of the collector surface
$S$ absorbed solar radiation, $S \approx 0.96 \tau \alpha_{beam} I_T$ also from Penn State, where $I_T$ is the incident solar radiation, I could just assume a panel's absorption efficiency rating
$U_L$ is the total loss coefficient through the panels etc (what about to the fluid)
$T_i$ is the fluid inlet temperature
$T_a$ is the ambient air temperature

Assuming glass panel
Chp 6 in book above
eqn 6.4.9 gives an empirical equation for the top loss coefficient of a solar collector
6.4.10 gives bottom loss coefficient equation, edge loss can be neglected, total loss is just summed
6.9.4 and 6.4.9 can be solved iteratively with a guessed mean plate temperature (generally 10C higher than inlet fluid temperature)
6.5.15 and 6.5.16 give equations for plate to fluid

**Approach:**
I think I can correctly calculate the loss coefficients without having to assume too many things (but some)

# CV2 (Pump) 

# CV3 (Storage Tank)
Called a thermal storage tank, storage is more efficient when stratified and *not* mixed [collection of papers](https://www.sciencedirect.com/topics/engineering/tank-thermal-energy-storage) 