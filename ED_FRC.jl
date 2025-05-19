# Packages
using JuMP
using LinearAlgebra
using MosekTools # Probar otro con SOC (Gurobi)
using Gurobi
using SCS

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
# N           = [   1,   30,   30,   25,   10,   12,   24,   16,    4,   19,   10,   21,   11,   13,   14,    8,    3,    7,   13,   16,    4,   19,   20,   21,   11,   13,   14,    8,    3,    7] #          - Number of units
# Pg_cost     = [  10,   95,   50,   45,   30,   75,   60,   80,   82,   58,   70,   45,   30,  175,  100,   80,   22,   38,   60,   80,   82,   58,   70,   45,   30,  175,  100,   80,   22,   38] # €/MWh    - Marginal cost
# Pg_nl_cost  = [   0,  500,  500,  250,  290,  380,  150,  300,  250,  285,  140,  250,  290,  380,  250,  200,  425,  385,  450,  620,  100,  185,  300,  260,  290,  180,  350,  220,  425,  385] # €        - No-load cost
# Pg_Gen_lb   = [1800,  250,   75,  100,  120,  175,  210,  110,   30,  185,   75,   80,  120,  175,  100,  180,  100,  125,  130,  110,  120,  185,   85,   80,  120,  195,  110,   80,  260,  125] # MW       - Power lower bound
# Pg_Gen_ub   = [1800,  500,  150,  400,  500,  550,  600,  450,  350,  740,  150,  400,  500,  650,  400,  750,  550,  330,  600,  450,  350,  740,  150,  400,  500,  350,  400,  350,  450,  240] # MW       - Power upper bound
# Rg_max      = [   0,  185,   50,   75,   65,   70,  120,   60,   90,  120,   50,   75,   65,   70,  120,   60,   90,  120,  120,   60,   90,  120,   50,   75,   65,   70,  120,   60,   90,  120] # MW       - FR provision
# Hg          = [   5,    5,    5,    5,    6,    7,    5,    6,    5,    6,    5,    5,    6,    7,    5,    6,    4,    6,    5,    6,    5,    6,    5,    8,    6,    7,    5,    6,    4,    6] # s        - Inertia constant

N           = [   1,   30,   30] #          - Number of units
Pg_cost     = [  10,   95,   50] # €/MWh    - Marginal cost
Pg_nl_cost  = [   0,  450,  500] # €        - No-load cost
Pg_Gen_lb   = [1800,  250,   75] # MW       - Power lower bound
Pg_Gen_ub   = [1800,  500,  150] # MW       - Power upper bound
Rg_max      = [   0,  185,   50] # MW       - FR provision
Hg          = [   5,    5,    5] # s        - Inertia constant

Tg          = 8                  # s        - FR delivery time

# Value of Lost Load (VoLL)
VoLL = 25000 # €/MWh

# Frecuency data
f_0         = 50 # Hz
Δf_max      = 0.8 # Hz
Δfss_max    = 0.5 # Hz

# Power from Renewable Energy Source (RES)
P_RES       = 40 *10^3 # MW (Max power installed)
# Capacity factor of RES (<= 1)
cf_RES      = [25, 25, 25, 25] ./ 100
Rs_ub       = 0 # MW (EFR)
Ts          = 1 # s

# Maximum power infeed
P_L_max = maximum(Pg_Gen_ub)
H_l = Hg[argmax(Pg_Gen_ub)]

# Load damping
D = 1.5 *10^-2 # %/Hz

# Total of power demand
Pd = [15, 30] .*10^3 # MW

# Number of periods
T = length(Pd)

# Check if the length of all data are the same
if (length(Pg_cost) == length(Pg_nl_cost) == length(Pg_Gen_lb) == length(Pg_Gen_ub) == length(Rg_max) == length(Hg))
    gTypes = length(Pg_cost)
else
    println("ERROR: Wrong data input\n")
end


########## Model creation ##########
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "OutputFlag", 0)


########## Variables ##########
@variable(model, Ng[1:gTypes, 1:T], Int) # Number of generators of each gTypes
@variable(model, Pg[1:gTypes, 1:T]) # Power from each generator
@variable(model, Rg[1:gTypes, 1:T]) # PRF from each generators
@variable(model, P_curt[1:T]) # RES Power curtailment
@variable(model, Rs[1:T]) # EFR from BESS (Battery Energy Storage Systems)
@variable(model, H[1:T]) # System inertia
@variable(model, P_L[1:T]) # Power infeed
@variable(model, Dshed[1:T] >= 0) # Demand shedding


########## Objective ##########
@objective(model, Min, sum(Ng[g, t] * Pg_nl_cost[g] + Pg[g, t] * Pg_cost[g] for g in 1:gTypes, t in 1:T) + sum(VoLL * Dshed[t] for t in 1:T))

########## Constraints ##########
@constraint(model, [t in 1:T], sum(Pg[g, t] for g in 1:gTypes) + P_RES * cf_RES[t] - P_curt[t] == Pd[t] - Dshed[t])

# Number of active generators of each group constraint
@constraint(model, [g in 1:gTypes, t in 1:T], Ng[g, t] >= 0)
@constraint(model, [g in 1:gTypes, t in 1:T], Ng[g, t] <= N[g])
@constraint(model, [t in 1:T], Ng[1, t] == N[1])

# Power from each generator constraints
@constraint(model, [g in 1:gTypes, t in 1:T], Pg[g, t] >= Pg_Gen_lb[g] * Ng[g, t])
@constraint(model, [g in 1:gTypes, t in 1:T], Pg[g, t] <= Pg_Gen_ub[g] * Ng[g, t])

# PRF provision from g constraints
@constraint(model, [g in 1:gTypes, t in 1:T], Rg[g, t] >= 0)
@constraint(model, [g in 1:gTypes, t in 1:T], Rg[g, t] <= Rg_max[g] * Ng[g, t])
@constraint(model, [g in 1:gTypes, t in 1:T], Rg[g, t] <= Pg_Gen_ub[g] * Ng[g, t] - Pg[g, t])

# EFR provision from RES constraints
@constraint(model, [t in 1:T], P_curt[t] >= 0)
@constraint(model, [t in 1:T], P_curt[t] <= P_RES * cf_RES[t])
@constraint(model, [t in 1:T], Rs[t] >= 0) # Rs (storage)
@constraint(model, [t in 1:T], Rs[t] <= Rs_ub)


# Largest power infeed constraints
@constraint(model, [t in 1:T], P_L[t] <= P_L_max)
@constraint(model, [g in 1:gTypes, t in 1:T], Pg[g, t] <= P_L[t] * Ng[g, t])
# McCormick relaxation
# @variable(model, aux[1:gTypes])
# @constraint(model, [i in 1:gTypes], aux[i] >= 0)                                                # (0 * Ng[i] + P_L * 0 - 0 * 0)
# @constraint(model, [i in 1:gTypes], aux[i] >= P_L_max * Ng[i] + P_L * N[i] - P_L_max * N[i])
# @constraint(model, [i in 1:gTypes], aux[i] <= P_L * N[i])                                       # (0 * Ng[i] + P_L * N[i] - 0 * N[i])
# @constraint(model, [i in 1:gTypes], aux[i] <= P_L_max * Ng[i])                                  # (P_L_max * Ng[i] + P_L * 0 - P_L_max * 0)
# @constraint(model, [i in 1:gTypes], Pg[i] <= aux[i])


# Equation 8: System inertia
@constraint(model, [t in 1:T], H[t] == sum(Hg[g] * Pg_Gen_ub[g] * Ng[g, t] for g in 1:gTypes) - P_L[t] * H_l)

# Equation 9: Quasi-steady-state security constraint
@constraint(model, [t in 1:T], Rs[t] + sum(Rg[g, t] for g in 1:gTypes) >= P_L[t] - Δfss_max * D * Pd[t])

# Equiation 13: Nadir constraint without load damping
# @constraint(model, (H / f_0 - Rs * Ts / (4*Δf_max)) * sum(Rg) >= (P_L - Rs)^2 * Tg / (4*Δf_max))

# Equation 19: Nadir constraint with load damping
# @constraint(model, (H / f_0 - Rs * Ts / (4*Δf_max)) * sum(Rg) >= (P_L - Rs)^2 * Tg / (4*Δf_max) - (P_L - Rs) * Tg * D * Pd / 4)




# RotatedSecondOrderCone()
# 2x₁x₂ ≥ ∣∣x₃₋ₙ∣∣²
# @constraint(model, [x₁, x₂, x₃, ..., xₙ] in RotatedSecondOrderCone())


# Auxiliar variables
@variable(model, x[1:5, 1:T])
@constraint(model, [t in 1:T], x[1, t] == H[t]/f_0 - Rs[t] * Ts / (4*Δf_max))
@constraint(model, [t in 1:T], x[2, t] == sum(Rg[g, t] for g in 1:gTypes) / Tg)
@constraint(model, [t in 1:T], x[3, t] == (P_L[t] - Rs[t]) / sqrt(4*Δf_max))
@constraint(model, [t in 1:T], x[4, t] == D * Pd[t])
@constraint(model, [t in 1:T], x[5, t] == (P_L[t] - Rs[t]) / 4)


# Equiation 13: Nadir constraint without load damping using RotatedSecondOrderCone()
# x1*x2 >= x3^2
# @constraint(model, [x[1,1]/2, x[2,1], x[3,1]] in RotatedSecondOrderCone())


# Equation 19: Nadir constraint with load damping using RotatedSecondOrderCone()
@variable(model, s_aux[1:T] >= 0)
@variable(model, t_aux[1:T] >= 0)
# x1*x2 >= x3^2 - x5*x4
# x1*x2 = s^2
@constraint(model, [t in 1:T], [x[1, t]/2, x[2, t], s_aux[t]] in RotatedSecondOrderCone())
# x5*x4 = t^2
@constraint(model, [t in 1:T], [x[5, t]/2, x[2, t], t_aux[t]] in RotatedSecondOrderCone())
# s^2 >= x3^2 - t^2 -----> s^2 + t^2 >= x3^2
@constraint(model, [t in 1:T], [s_aux[t], t_aux[t], x[3, t]] in RotatedSecondOrderCone())


clearTerminal()
########## Optimization ##########
optimize!(model)


########## Solution display ##########
# clearTerminal()
println("Capacidad = $(sum(Pg_Gen_ub[n] for n in 1:gTypes))")
if termination_status(model) == OPTIMAL || termination_status(model) == LOCALLY_SOLVED
    
    println("\n\n##### $(termination_status(model)) solution found #####")

    for t in 1:T
        println("\nDemand = ", Pd[t], " MW")

        println("Generated power = ", round(sum(value(Pg[g, t]) for g in 1:gTypes), digits = 3), " MW")

        println("RES power supply = ", round(value(P_RES * cf_RES[t] - P_curt[t]), digits = 3), " MW")
        println("RES accommodated = ", round(value(P_curt[t]), digits = 2), " MW")
        println("Demmand shedding = ", round(value(Dshed[t]), digits = 2), " MW")
        

        # idx = cumsum([1; N[1:end-1]])
        for i in 1:gTypes
            println("")
            println("Number of generators of type $i = $(value(Ng[i, t]))")
            println("Power supplied by gen. type $i  = $(round(value(Pg[i, t]), digits = 2)) MW")
            println("PFR provision from gen units $i = $(round(value(Rg[i, t]), digits = 2)) MW")
            println("Operation cost = $(round(value(Pg[i, t]) * Pg_cost[i] / 1000, digits = 2)) k€")
        end
        println("")
        println("Load infeed = ", round(value(P_L[t]), digits = 3), " MW")
        println("")
        println("Global PRF provision = $(round(sum(value(Rg[i, t]) for i in 1:gTypes), digits = 2)) MW")
        cost_t = sum(value(Ng[g, t]) * Pg_nl_cost[g] + value(Pg[g, t]) * Pg_cost[g] for g in 1:gTypes) + VoLL * value(Dshed[t])
        println("Total cost period $t = $(round(cost_t/1000, digits = 3)) k€")
    end
else
    println("ERROR: ", termination_status(model))

end