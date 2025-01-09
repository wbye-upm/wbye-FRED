# Packages
using JuMP
using Gurobi # Probar otro con SOC (Gurobi)

# Function to clean the terminal
function clearTerminal()

    # For Windows terminal
    if Sys.iswindows()
        Base.run(`cmd /c cls`)

    # Other terminals based on Unix
    else
        Base.run(`clear`)
    end
    
end

########## Data ##########
# Generators data
N           = [   1,   30,   30] #          - Number of units
Pg_cost     = [  10,   95,   50] # €/MWh    - Marginal cost
Pg_nl_cost  = [   0,  500,  500] # €        - No-load cost
Pg_Gen_lb   = [1800,  250,   75] # MW       - Power lower bound
Pg_Gen_ub   = [1800,  500,  150] # MW       - Power upper bound
Rg_max      = [   0,  185,   50] # MW       - FR provision
Hg          = [   5,    5,    5] # s        - Inertia constant
Tg          = 8                  # s        - FR delivery time

# Frecuency data
f_0         = 50 # Hz
Δf_max      = 0.8 # Hz
Δfss_max    = 0.5 # Hz

# Power from Renewable Energy Source (RES)
P_RES       = 40 *10^3 # MW (Max power installed)
cf_RES      = 0.25 # Capacity factor of RES (<= 1)
Rs_ub       = 0 # MW (EFR)
Ts          = 1 # s

# Maximum power infeed
P_L_max = maximum(Pg_Gen_ub)
H_l = Hg[argmax(Pg_Gen_ub)]

# Load damping
D = 1.5 # %/Hz

# Total of power demand
Pd = 24 *10^3 # MW


# Number of total generators
# G = sum(N)


# Check if the length of all data are the same
if (length(Pg_cost) == length(Pg_nl_cost) == length(Pg_Gen_lb) == length(Pg_Gen_ub) == length(Rg_max) == length(Hg))
    gTypes = length(Pg_cost)
else
    println("ERROR: Wrong data input")
end


########## Model creation ##########
model = Model(Gurobi.Optimizer)
# set_silent(model)


########## Variables ##########
@variable(model, Ng[1:gTypes], Int) # Number of generators of each gTypes
@variable(model, Pg[1:gTypes]) # Power from each generator
@variable(model, Rg[1:gTypes]) # PRF from each generators
@variable(model, P_curt) # RES Power curtailment
@variable(model, Rs) # EFR from BESS (Battery Energy Storage Systems)
@variable(model, H) # System inertia
@variable(model, P_L) # Power infeed


########## Objective ##########
# @objective(model, Min, sum(Pg[n] * Pg_cost[n] for n in 1:G) + VoLL * loadShedding) 
    # VoLL = Value of Loss Load (Valor alto, ej. 20000 €/MWh)
    # loadShedding = Power demand reduction for fullfilling the lack of power generation with high cost
# Objetive without considering load shedding
# @objective(model, Min, sum(Pg[n] * Pg_cost[n] for n in 1:G))
@objective(model, Min, sum(Ng[n] * Pg_nl_cost[n] + Pg[n] * Pg_cost[n] for n in 1:gTypes))


########## Constraints ##########
# @constraint(model, sum(Pg) + P_RES - P_curt == Pd - loadShedding)
@constraint(model, sum(Pg) + P_RES * cf_RES - P_curt == Pd)

# # Number of active generators of each group constraint
@constraint(model, [i in 1:gTypes], 0 <= Ng[i])
@constraint(model, [i in 1:gTypes], Ng[i] <= N[i])
@constraint(model, Ng[1] == N[1])

# Power from each generator constraints
@constraint(model, [i in 1:gTypes], Pg_Gen_lb[i] * Ng[i] <= Pg[i])
@constraint(model, [i in 1:gTypes], Pg[i] <= Pg_Gen_ub[i] * Ng[i])

# PRF provision from g constraints
@constraint(model, [i in 1:gTypes], 0 <= Rg[i])
@constraint(model, [i in 1:gTypes], Rg[i] <= Rg_max[i] * Ng[i])
@constraint(model, [i in 1:gTypes], Rg[i] <= Pg_Gen_ub[i] * Ng[i] - Pg[i])

# EFR provision from RES constraints
@constraint(model, 0 <= P_curt)
@constraint(model, P_curt <= P_RES * cf_RES)
@constraint(model, 0 <= Rs) # Rs (storage)
@constraint(model, Rs <= Rs_ub)


# Largest power infeed constraints
@constraint(model, P_L <= P_L_max)
# @constraint(model, [i in 1:gTypes], Pg[i] <= P_L * Ng[i])
@variable(model, aux[1:gTypes])
@constraint(model, [i in 1:gTypes], aux[i] >= 0 * Ng[i] + P_L * 0 - 0 * 0)
@constraint(model, [i in 1:gTypes], aux[i] >= P_L_max * Ng[i] + P_L * N[i] - P_L_max * N[i])
@constraint(model, [i in 1:gTypes], aux[i] <= 0 * Ng[i] + P_L * N[i] - 0 * N[i])
@constraint(model, [i in 1:gTypes], aux[i] <= P_L_max * Ng[i] + P_L * 0 - P_L_max * 0)
@constraint(model, [i in 1:gTypes], Pg[i] <= aux[i])



# Equation 8: System inertia
@constraint(model, H == sum(Hg .* Pg_Gen_ub .* Ng) - P_L * H_l)

# Equation 9: Quasi-steady-state security constraint
# @constraint(model, (P_L - Rs - sum(Rg) / (D * Pd)) <= Δfss_max)
@constraint(model, Rs + sum(Rg) >= P_L - Δfss_max * D * Pd)

# Equiation 13: Nadir constraint without load damping
@constraint(model, (H / f_0 - Rs * Ts / (4*Δf_max)) * sum(Rg) >= (P_L - Rs)^2 * Tg / (4*Δf_max))
##revisar: Implementar un cono, llegar a "continuous convex" conic equation
##revisar: al meter como cono y con valores binarios tiene que salir MISOCP

# Equation 19: Nadir constraint with load damping
# @constraint(model, (H / f_0 - Rs * Ts / (4*Δf_max)) * sum(Rg) >= (P_L - Rs)^2 * Tg / (4*Δf_max) - (P_L - Rs) * Tg * D * Pd / 4)
##revisar: si hay algun caso que gurobi no pueda resolver no convexidades y haya que acudir a Ipopt (con relajaciones de las variables enteras)

##revisar: Boyd
clearTerminal()
########## Optimization ##########
optimize!(model)


########## Solution display ##########
# clearTerminal()

if termination_status(model) == OPTIMAL || termination_status(model) == LOCALLY_SOLVED
    println("\n\n##### $(termination_status(model)) solution found #####\n")
    println("H / f_0 - Rs * Ts / (4*Δf_max) = ", round(value((Hg' * Pg_Gen_ub) / f_0 - Rs * Ts / (4*Δf_max)), digits = 3))

    # for i in 1:G
    #     println("Rg del generador ", i, " = ", round(value(Rg[i]),digits = 3))
    # end
    println("RG = ", round(value(sum(Rg)), digits = 3))
    println("Rs = ", round(value(Rs), digits = 3))
    println("P_L - Δfss_max * D * Pd = ", round(value(P_L - Δfss_max * D * Pd), digits = 3))

    println("\nDemand = ", Pd, " MW")

    println("Generated power = ", round(sum(value(Pg[i]) for i in 1:gTypes), digits = 3), " MW")

    println("RES power supply = ", round(value(P_RES * cf_RES - P_curt), digits = 3), " MW")
    println("RES accommodated = ", round(value(P_curt), digits = 3), " MW")

    # for i in 1:G
    #     println("Potencia generador ", i, " : ", round(value(Pg[i]), digits = 3), " MW")
    # end

    # print("\n", Hg, "\n\n")
    

    idx = cumsum([1; N[1:end-1]])
    for i in 1:gTypes
        println("")
        println("Number of generators of type $i = $(Int(value(Ng[i])))")
        println("Power supplied by gen. type $i  = $(round(value(Pg[i]), digits = 3)) MW")
        println("PRF provision from gen units $i = $(round(value(Rg[i]), digits = 3)) MW")
        println("Operation cost = $(round(value(Pg[i]) * Pg_cost[i] / 1000, digits = 3)) k€")
    end
    println("")
    println("Load infeed = ", round(value(P_L), digits = 3), " MW")
    println("")
    println("Total cost = $(round(objective_value(model), digits = 3)/1000) k€")

else
    println("ERROR: ", termination_status(model))

end