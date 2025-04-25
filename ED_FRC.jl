# Packages
using JuMP
using LinearAlgebra
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

N           = [   1,   30,   30,   50] #          - Number of units
Pg_cost     = [  10,   95,   50,   50] # €/MWh    - Marginal cost
Pg_nl_cost  = [   0,  500,  500,  550] # €        - No-load cost
Pg_st_cost  = [   0,   10,   20,   15] # €        - Startup cost
Pg_Gen_lb   = [1800,  250,  100,  150] # MW       - Power lower bound
Pg_Gen_ub   = [1800,  500,  300,  450] # MW       - Power upper bound
Rg_max      = [   0,  185,   50,  100] # MW       - FR provision
Hg          = [   5,    5,    5,    5] # s        - Inertia constant
Pgrr        = [   0,   50,   50,   50] # MW/h     - Maximum ramp rate of each generator
Tst         = [ 999,    2,    3,    1] # h        - Generators startup time
Tmdt        = [ 999,    1,    2,    1] # h        - Generators minimum down time
Tmut        = [ 999,    3,    1,    2] # h        - Generators minimum up time

Tg          = 8                        # s        - FR delivery time

Ng_init     = [   1,   20,   16,  50] #           - Number of starting units
Pg_init     = [1800, 5000, 4000, 9000] #          - Starting power of each cluster

# Value of Lost Load (VoLL)
VoLL = 25000 # €/MWh

# Frecuency data
f_0         = 50 # Hz
Δf_max      = 0.8 # Hz
Δfss_max    = 0.5 # Hz

# Total of power demand in MW
# Pd          = [24, 22, 20, 17, 15, 17, 21, 22, 23, 25, 27, 26, 30, 27, 25] .* 1000
Pd          = [25, 20, 27] .* 1000

# Power from Renewable Energy Source (RES)
P_RES       = 40 *10^3 # MW (Max power installed)
# Capacity factor of RES (<= 1)
# cf_RES      = [20, 17, 26, 27, 30, 18, 16, 10,  9, 25, 22, 26, 28, 30, 15] ./ 100
cf_RES      = [00, 00, 00] ./ 100
Rs_ub       = 0 # MW (EFR)
Ts          = 1 # s

# Maximum power infeed
P_L_max = maximum(Pg_Gen_ub)
H_l = Hg[argmax(Pg_Gen_ub)]

# Load damping
D = 1.5 *10^-2 # %/Hz


T = length(Pd)


# Check if the length of all data are the same
if (length(Pg_cost) == length(Pg_nl_cost) == length(Pg_Gen_lb) == length(Pg_Gen_ub) == length(Rg_max) == length(Hg))
    gTypes = length(Pg_cost)
else
    println("ERROR: Wrong data input\n")
end


########## Model creation ##########
model = Model(Gurobi.Optimizer)
# if string(typeof(backend(model).optimizer.model.optimizer)) == "Gurobi.Optimizer"
#     set_optimizer_attribute(model, "OutputFlag", 0)
# elseif string(typeof(backend(model).optimizer.model.optimizer)) == "SCS.Optimizer"
#     set_optimizer_attribute(model, "verbose", 0)
# end
set_optimizer_attribute(model, "OutputFlag", 0)


########## Variables ##########
@variable(model, Ng[1:gTypes], Int) # Number of generators of each gTypes
@variable(model, Pg[1:gTypes]) # Power from each generator
@variable(model, Rg[1:gTypes]) # PRF from each generators
@variable(model, P_curt) # RES Power curtailment
@variable(model, Rs) # EFR from BESS (Battery Energy Storage Systems)
@variable(model, H) # System inertia
@variable(model, P_L) # Power infeed
@variable(model, Dshed >= 0) # Demand shedding

@variable(model, cf_RES_var) # Auxiliar variable for cf_RES
@variable(model, Pd_var) # Auxiliar variable for Pd
@variable(model, Ng_var[1:gTypes]) # Auxiliar variable for Ng[t-1]
@variable(model, Pg_var[1:gTypes]) # Auxiliar variable for Pg[t-1]


########## Objective ##########
@objective(model, Min, sum(Ng[n] * Pg_nl_cost[n] + Pg[n] * Pg_cost[n] for n in 1:gTypes) + VoLL * Dshed)


########## Constraints ##########
@constraint(model, sum(Pg) + P_RES * cf_RES_var - P_curt >= Pd_var - Dshed)

# Number of active generators of each group constraint
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
@constraint(model, P_curt <= P_RES * cf_RES_var)
@constraint(model, 0 <= Rs) # Rs (storage)
@constraint(model, Rs <= Rs_ub)


# Largest power infeed constraints
@constraint(model, P_L <= P_L_max)

# @constraint(model, [i in 1:gTypes], Pg[i] <= P_L * Ng[i])
# McCormick relaxation
@variable(model, aux[1:gTypes])
@constraint(model, [i in 1:gTypes], aux[i] >= 0)                                                # (0 * Ng[i] + P_L * 0 - 0 * 0)
@constraint(model, [i in 1:gTypes], aux[i] >= P_L_max * Ng[i] + P_L * N[i] - P_L_max * N[i])
@constraint(model, [i in 1:gTypes], aux[i] <= P_L * N[i])                                       # (0 * Ng[i] + P_L * N[i] - 0 * N[i])
@constraint(model, [i in 1:gTypes], aux[i] <= P_L_max * Ng[i])                                  # (P_L_max * Ng[i] + P_L * 0 - P_L_max * 0)
@constraint(model, [i in 1:gTypes], Pg[i] <= aux[i])


# Equation 8: System inertia
@constraint(model, H == sum(Hg .* Pg_Gen_ub .* Ng) - P_L * H_l)

# Equation 9: Quasi-steady-state security constraint
@constraint(model, Rs + sum(Rg) >= P_L - Δfss_max * D * Pd_var)

# Equiation 13: Nadir constraint without load damping
# @constraint(model, (H / f_0 - Rs * Ts / (4*Δf_max)) * sum(Rg) >= (P_L - Rs)^2 * Tg / (4*Δf_max))

# Equation 19: Nadir constraint with load damping
# @constraint(model, (H / f_0 - Rs * Ts / (4*Δf_max)) * sum(Rg) >= (P_L - Rs)^2 * Tg / (4*Δf_max) - (P_L - Rs) * Tg * D * Pd_var / 4)



# RotatedSecondOrderCone()
# 2x₁x₂ ≥ ∣∣x₃₋ₙ∣∣²
# @constraint(model, [x₁, x₂, x₃, ..., xₙ] in RotatedSecondOrderCone())


# Auxiliar variables
@variable(model, x[1:5, 1])
@constraint(model, x[1,1] == H/f_0 - Rs * Ts / (4*Δf_max))
@constraint(model, x[2,1] == sum(Rg) / Tg)
@constraint(model, x[3,1] == (P_L - Rs) / sqrt(4*Δf_max))
@constraint(model, x[4,1] == D * Pd_var)
@constraint(model, x[5,1] == (P_L - Rs) / 4)


# Equiation 13: Nadir constraint without load damping using RotatedSecondOrderCone()
# x1*x2 >= x3^2
# @constraint(model, [x[1,1]/2, x[2,1], x[3,1]] in RotatedSecondOrderCone())


# Equation 19: Nadir constraint with load damping using RotatedSecondOrderCone()
@variable(model, s >= 0)
@variable(model, t >= 0)
# x1*x2 >= x3^2 - x5*x4
# x1*x2 = s^2
@constraint(model, [x[1,1]/2, x[2,1], s] in RotatedSecondOrderCone())
# x5*x4 = t^2
@constraint(model, [x[5,1]/2, x[2,1], t] in RotatedSecondOrderCone())
# s^2 >= x3^2 - t^2 -----> s^2 + t^2 >= x3^2
@constraint(model, [s, t, x[3,1]] in RotatedSecondOrderCone())


# Ramp limits
# Max variation of the power generated by each generator
# @constraint(model, [i in 1:gTypes], -Pgrr[i] * Ng_buffer[nPeriod - 1][i] <= Pg[i] - Pg_buffer[nPeriod - 1][i])
@constraint(model, [i in 1:gTypes], -Pgrr[i] * Ng_var[i] <= Pg[i] - Pg_var[i])
# # @constraint(model, [i in 1:gTypes], Pg[i] - Pg_buffer[nPeriod - 1][i] <= Pgrr[i] * Ng_buffer[nPeriod][i])
@constraint(model, [i in 1:gTypes], Pg[i] - Pg_var[i] <= Pgrr[i] * Ng[i])



clearTerminal()

########## Period dependant constraints update ##########
Ng_buffer = []
Pg_buffer = []
for nPeriod in 1:lastindex(Pd)

    fix(cf_RES_var, cf_RES[nPeriod]; force = true)
    fix(Pd_var, Pd[nPeriod]; force = true)

    if nPeriod == 1
        fix.(Ng_var, Ng_init; force = true)
        fix.(Pg_var, Pg_init; force = true)
    else
        fix.(Ng_var, Ng_buffer[nPeriod-1]; force = true)
        fix.(Pg_var, Pg_buffer[nPeriod-1]; force = true)
    end

    ########## Optimization ##########
    optimize!(model)

    println("\n")
    println("Period: ", nPeriod)
    for i in 1:gTypes
        println("")
        println("Number of generators of type $i = $((value(Ng[i])))")
        println("Power supplied by gen. type $i  = $(round(value(Pg[i]), digits = 2)) MW")
        println("PFR provision from gen units $i = $(round(value(Rg[i]), digits = 2)) MW")
        println("Operation cost = $(round(value(Pg[i]) * Pg_cost[i] / 1000, digits = 2)) k€")
    end
    push!(Ng_buffer, round.(value.(Ng), digits = 0))
    push!(Pg_buffer, round.(value.(Pg), digits = 2))
    println("\nGenerated power = ", round(sum(value(Pg[i]) for i in 1:gTypes), digits = 3), " MW")
    println("Demand shedding = ", round(value(Dshed), digits = 3), " MW")
    println("Total cost = $(round(objective_value(model)/1000, digits = 3)) k€")

end


########## Solution display ##########
# clearTerminal()
# println("Capacidad = $(sum(Pg_Gen_ub[n] for n in 1:gTypes))")
# if termination_status(model) == OPTIMAL || termination_status(model) == LOCALLY_SOLVED
#     println("\n\n##### $(termination_status(model)) solution found #####\n")
#     println("H / f_0 - Rs * Ts / (4*Δf_max) = ", round(value((Hg' * Pg_Gen_ub) / f_0 - Rs * Ts / (4*Δf_max)), digits = 3))

#     # for i in 1:G
#     #     println("Rg del generador ", i, " = ", round(value(Rg[i]),digits = 3))
#     # end
#     println("RG = ", round(value(sum(Rg)), digits = 3))
#     println("Rs = ", round(value(Rs), digits = 3))
#     println("P_L - Δfss_max * D * Pd = ", round(value(P_L - Δfss_max * D * Pd[1]), digits = 3))

#     println("\nDemand = ", Pd[1], " MW")

#     println("Generated power = ", round(sum(value(Pg[i]) for i in 1:gTypes), digits = 3), " MW")

#     println("RES power supply = ", round(value(P_RES * cf_RES[1] - P_curt), digits = 3), " MW")
#     println("RES accommodated = ", round(value(P_curt), digits = 2), " MW")
    

#     # idx = cumsum([1; N[1:end-1]])
#     for i in 1:gTypes
#         println("")
#         println("Number of generators of type $i = $(round(value(Ng[i]), digits = 0))")
#         println("Power supplied by gen. type $i  = $(round(value(Pg[i]), digits = 2)) MW")
#         println("PFR provision from gen units $i = $(round(value(Rg[i]), digits = 2)) MW")
#         println("Operation cost = $(round(value(Pg[i]) * Pg_cost[i] / 1000, digits = 2)) k€")
#     end
#     println("")
#     println("Load infeed = ", round(value(P_L), digits = 3), " MW")
#     println("")
#     println("Demand shedding = ", round(value(Dshed), digits = 3), " MW")
#     println("Global PRF provision = $(round(sum(value(Rg[i]) for i in 1:gTypes), digits = 2)) MW")
#     println("Total cost = $(round(objective_value(model)/1000, digits = 3)) k€")

# else
#     println("ERROR: ", termination_status(model))

# end