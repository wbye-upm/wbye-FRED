# Packages
using JuMP
using LinearAlgebra
using Gurobi

using DataFrames
using Plots

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

# https://www.esios.ree.es/es/unidades-fisicas
########## Data ##########
Types       = [ "Nuclear", "CC_1", "CC_2", "CC_3", "CC_4", "Carbón_1", "Carbón_2", "COGEN_1", "COGEN_2", "COGEN_3", "COGEN_4", "COGEN_5", "Fuel", "Derivados"]
N           = [         7,      5,     33,     18,     20,          2,          1,        97,        58,        86,        31,        12,      2,          52] #        - Number of units
C_m         = [        11,   47.5,   42.5,     55,   52.5,         39,         34,     107.5,      92.5,        85,      81.5,        70,    105,        97.5] # €/MWh  - Marginal cost
C_nl        = [      22.5,      9,   12.5,     16,     19,         38,         45,       6.5,       8.5,        10,      12.5,        17,      9,          12] # €/h    - No-load cost
Pg_Gen_lb   = [     457.5,   25.2,  120.2,  174.5,    187,        154,        256,      0.51,      1.61,      4.13,     10.32,     41.99,    1.8,         6.3] # MW     - Power lower bound
Pg_Gen_ub   = [    1016.7,   55.9,  267.1,  387.7,  415.3,     343.95,        570,      1.14,      3.57,      9.18,     22.94,      91.1,   3.95,       13.96] # MW     - Power upper bound
Pg_rr       = [      14.5,   22.5,   42.5,     60,     85,         75,      112.5,      12.5,        19,      28.5,      37.5,        50,     12,          21] # MW/h   - Maximum ramp rate for each generator of the cluster
Hg          = [         9,      6,      6,      7,      7,          8,          8,         3,         3,         4,         5,         5,      5,           5] # s      - Inertia constant
T_st        = [        48,      2,      2,      2,      2,         12,         15,         0,         0,         0,         1,         2,      1,           2] # h      - Startup time
T_mut       = [        48,      4,      4,      4,      4,          6,          6,         2,         3,         3,         3,         4,      1,           2] # h      - Minimum up time
T_mdt       = [        48,      4,      4,      4,      4,          6,          6,         2,         2,         2,         2,         2,      1,           1] # h      - Minimum down time

Rg_max      = Pg_Gen_ub.*0.05   # MW    - FR provision
Tg          = 8                 # s     - FR delivery time

# https://demanda.ree.es/visiona/peninsula/demandaau/acumulada/2025-01-20
# Hora:         00:00   01:00   02:00   03:00   04:00   05:00   06:00   07:00   08:00   09:00   10:00   11:00   12:00   13:00   14:00   15:00   16:00   17:00   18:00   19:00   20:00   21:00   22:00   23:00
# Total of power demand
Pd = [27.062, 24.434, 22.861, 22.076, 21.858, 22.104, 24.402, 29.278, 34.13, 36.113, 35.758, 35.066, 34.56, 34.924, 35.315, 35.564, 36.086, 36.699, 37.582, 38.378, 38.714, 38.521, 35.598, 31.529] .*10^3 # MW
# Capacity factor of RES
cf_RES      = [12.264, 11.9488, 12.1419, 12.8872, 12.636, 12.3849, 12.6209, 12.2465, 12.0326, 12.5942, 14.1442, 14.9837, 14.8953, 15.3384, 14.4291, 13.0651, 11.4547, 10.093, 10.7326, 11.8674, 11.85, 11.6837, 11.264, 11.05] ./ 100

# Power from Renewable Energy Source (RES)
P_RES       = 86 *10^3 # MW (Max power installed)

Rs_ub       = 0 # MW (EFR)
Ts          = 1 # s

# Max Rate of Change of Frecuency (RoCoF)
RoCoF_max = 1 # Hz/s

# Number of periods
T = length(Pd)

# Value of Lost Load (VoLL)
VoLL = 250000 # €/MWh

# Frecuency data
f_0         = 50 # Hz
Δf_max      = 0.8 # Hz
Δfss_max    = 0.2 # Hz

# Maximum power infeed
P_L_max = maximum(Pg_Gen_ub)
H_l = Hg[argmax(Pg_Gen_ub)]

# Load damping
D = 1.5 *10^-2 # %/Hz


# Check if the length of all data are the same
gTypes = length(Types)


########## Model creation ##########
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "OutputFlag", 0) # Se silencian las salidas pro defecto del solver Gurobi


########## Variables ##########
@variable(model, Ng[1:gTypes, 1:T] >= 0, Int) # Number of generators of each gTypes
@variable(model, Ng_sd[1:gTypes, 1:T] >= 0, Int) # Number of generators of each gTypes - Shutdown
@variable(model, Ng_sg[1:gTypes, 1:T] >= 0, Int) # Number of generators of each gTypes - Start generating
@variable(model, Ng_st[1:gTypes, 1:T] >= 0, Int) # Number of generators set to Startup
@variable(model, Pg[1:gTypes, 1:T] >= 0) # Power from each cluster of generators
@variable(model, Rg[1:gTypes, 1:T] >= 0) # PRF from each generators
@variable(model, P_curt[1:T] >= 0) # RES Power curtailment
@variable(model, Rs[1:T] >= 0) # EFR from BESS (Battery Energy Storage Systems)
@variable(model, H[1:T] >= 0) # System inertia
@variable(model, P_L[1:T] >= 0) # Power infeed
@variable(model, Dshed[1:T] >= 0) # Demand shedding


########## Objective ##########
@objective(model, Min, sum(Ng[g, t] * C_nl[g] + Pg[g, t] * C_m[g] for g in 1:gTypes, t in 1:T) + sum(VoLL * Dshed[t] for t in 1:T))

########## Constraints ##########
@constraint(model, [t in 1:T], sum(Pg[g, t] for g in 1:gTypes) + P_RES * cf_RES[t] - P_curt[t] == Pd[t] - Dshed[t])

# Number of active generators of each group constraint
@constraint(model, [g in 1:gTypes, t in 1:T], Ng[g, t] <= N[g])
@constraint(model, [t in 1:T], Ng[1, t] == N[1])

# Power from each generator constraints
@constraint(model, [g in 1:gTypes, t in 1:T], Pg[g, t] >= Pg_Gen_lb[g] * Ng[g, t])
@constraint(model, [g in 1:gTypes, t in 1:T], Pg[g, t] <= Pg_Gen_ub[g] * Ng[g, t])

# Largest power infeed constraints
@constraint(model, [t in 1:T], P_L[t] <= P_L_max)
@constraint(model, [g in 1:gTypes, t in 1:T], Pg[g, t] <= P_L[t] * Ng[g, t])

# PFR provision from g constraints
@constraint(model, [g in 1:gTypes, t in 1:T], Rg[g, t] <= Rg_max[g] * Ng[g, t])
@constraint(model, [g in 1:gTypes, t in 1:T], Rg[g, t] <= Pg_Gen_ub[g] * Ng[g, t] - Pg[g, t])

# EFR provision from RES constraints
@constraint(model, [t in 1:T], P_curt[t] >= 0)
@constraint(model, [t in 1:T], P_curt[t] <= P_RES * cf_RES[t])
@constraint(model, [t in 1:T], Rs[t] >= 0) # Rs (storage)
@constraint(model, [t in 1:T], Rs[t] <= Rs_ub)

# RoCoF constraint
@constraint(model, [t in 1:T], P_L[t] * f_0 / (2*H[t]) <= RoCoF_max)

# System inertia
@constraint(model, [t in 1:T], H[t] == sum(Hg[g] * Pg_Gen_ub[g] * Ng[g, t] for g in 1:gTypes) - P_L[t] * H_l)

# Quasi-steady-state security constraint
@constraint(model, [t in 1:T], (Rs[t] + sum(Rg[g, t] for g in 1:gTypes) - P_L[t]) / (D * Pd[t]) <= Δfss_max)
@constraint(model, [t in 1:T], -(Rs[t] + sum(Rg[g, t] for g in 1:gTypes) - P_L[t]) / (D * Pd[t]) <= Δfss_max)


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

# Nadir constraint without load damping
# @constraint(model, (H / f_0 - Rs * Ts / (4*Δf_max)) * sum(Rg) >= (P_L - Rs)^2 * Tg / (4*Δf_max))

# Nadir constraint without load damping using RotatedSecondOrderCone()
# x1*x2 >= x3^2
# @constraint(model, [t in 1:T], [x[1,t]/2, x[2,t], x[3,t]] in RotatedSecondOrderCone())


# Nadir constraint with load damping
# @constraint(model, (H / f_0 - Rs * Ts / (4*Δf_max)) * sum(Rg) >= (P_L - Rs)^2 * Tg / (4*Δf_max) - (P_L - Rs) * Tg * D * Pd / 4)

# Nadir constraint with load damping using RotatedSecondOrderCone()
@variable(model, s_aux[1:T] >= 0)
@variable(model, t_aux[1:T] >= 0)
# x1*x2 >= x3^2 - x5*x4
# x1*x2 = s^2
@constraint(model, [t in 1:T], [x[1, t]/2, x[2, t], s_aux[t]] in RotatedSecondOrderCone())
# x5*x4 = t^2
@constraint(model, [t in 1:T], [x[5, t]/2, x[2, t], t_aux[t]] in RotatedSecondOrderCone())
# s^2 >= x3^2 - t^2 -----> s^2 + t^2 >= x3^2
@constraint(model, [t in 1:T], [s_aux[t], t_aux[t], x[3, t]] in RotatedSecondOrderCone())


for t in 1:T
    if t > 1
        # Ramp limits constraints
        @constraint(model, [g in 1:gTypes], -Pg_Gen_lb[g] * Ng_sd[g, t] - Pg_rr[g] * (Ng[g, t] - Ng_sg[g, t]) <= Pg[g, t] - Pg[g, t-1])
        @constraint(model, [g in 1:gTypes], Pg[g, t] - Pg[g, t-1] <= Pg_rr[g] * (Ng[g, t] - Ng_sg[g, t]) + Pg_Gen_lb[g] * Ng_sg[g, t])
        # Number of generators of each cluster variation
        @constraint(model, [g in 1:gTypes], Ng[g, t] == Ng[g, t-1] + Ng_sg[g, t] - Ng_sd[g, t])
    end
    for g in 1:gTypes
        # Start generating after startup
        if t > T_st[g]
            @constraint(model, Ng_sg[g, t] == Ng_st[g, t - T_st[g]])
        end
        # Minimum down time constraint
        if (t > T_mdt[g] && T_mdt[g] > 0) || (t > T_mdt[g]+1 && T_mdt[g] == 0)
            @constraint(model, Ng_st[g, t] <= N[g] - Ng[g, t-1] - sum(Ng_sd[g, i] for i in t-T_mdt[g]:t))
        elseif t-T_mdt[g] < 1
            @constraint(model, Ng_st[g, t] == 0)
        end
        # Minimun up time constraint
        if (t > T_mut[g] && T_mut[g] > 0) || (t > T_mut[g]+1 && T_mut[g] == 0)
            @constraint(model, Ng_sd[g, t] <= Ng[g, t-1] - sum(Ng_sg[g, i] for i in t-T_mut[g]:t-1))
        elseif t-T_mut[g] < 1
            @constraint(model, Ng_sd[g, t] == 0)
        end
    end
end


clearTerminal()
########## Optimization ##########
optimize!(model)


########## Solution display ##########
if termination_status(model) == OPTIMAL || termination_status(model) == LOCALLY_SOLVED
    
    println("\n\n##### $(termination_status(model)) solution found #####")

    for t in 1:T
        # println("\nDemand = ", Pd[t], " MW")

        println("Generated power = ", round(sum(value(Pg[g, t]) for g in 1:gTypes), digits = 3), " MW")

        println("RES power supply = ", round(value(P_RES * cf_RES[t] - P_curt[t]), digits = 3), " MW")
        println("RES accommodated = ", round(value(P_curt[t]), digits = 2), " MW")
        # println("Demmand shedding = ", round(value(Dshed[t]), digits = 2), " MW")
        println("\nPeriod $t:")
        for g in 1:gTypes
            # println("")
            println("Number of generators of type $g = $(value(Ng[g, t]))")
            println("Power supplied by gen. type $g  = $(round(value(Pg[g, t]), digits = 2)) MW")
            println("Power supplied by each gen per type $g  = $(round(value(Pg[g, t] / Ng[g, t]), digits = 2)) MW")
            # println("PFR provision from gen units $g = $(round(value(Rg[g, t]), digits = 2)) MW")
            # println("Operation cost = $(round(value(Pg[g, t]) * C_m[g] / 1000, digits = 2)) k€")
            # println("Ng_st = $(value(Ng_st[g, t]))")
            # println("Ng_sg = $(value(Ng_sg[g, t]))")
            # println("Ng_sd = $(value(Ng_sd[g, t]))")
        end
        # println("")
        println("Load infeed = ", round(value(P_L[t]), digits = 3), " MW")
        # println("")
        println("Global PRF provision = $(round(sum(value(Rg[g, t]) for g in 1:gTypes), digits = 2)) MW")
        cost_t = sum(value(Ng[g, t]) * C_nl[g] + value(Pg[g, t]) * C_m[g] for g in 1:gTypes) + VoLL * value(Dshed[t])
        # println("Total cost period $t = $(round(cost_t/1000, digits = 3)) k€")

        println("Δfss = ", value((Rs[t] + sum(Rg[g, t] for g in 1:gTypes) - P_L[t]) / (D * Pd[t])))
    end

    ########## Solution plotting ##########
    x = 1:T
    Pg_values = [value(Pg[g, t]/1000) for g in 1:gTypes, t in 1:T]
    y = Array{Float64}(undef, gTypes, T)
    y[1, :] .= Pg_values[1, :]
    plot(x, y[1, :], ylims=(0, maximum(Pd/1000)), fillrange=0, lw=0, label=Type[1], xlabel="Time [h]", ylabel="Power Output (MW)")
    for g in 2:gTypes
        y[g, :] .= y[g-1, :] .+ Pg_values[g, :]
        plot!(x, y[g, :], fillrange=y[g-1, :], lw=0, label=Type[g])
    end
    y_RES = [value(P_RES * cf_RES[t] - P_curt[t])/1000 for t in 1:T]
    plot!(x, y_RES + y[gTypes, :], fillrange=y[gTypes, :], lw=0, label="RES power")
    plot!(x, Pd/1000, label="Demand", color="red", lw=3)
    load_shedding = [value(Dshed[t]) for t in 1:T]
    display(current())

    println(solve_time(model))

    

else
    println("ERROR: ", termination_status(model))

end