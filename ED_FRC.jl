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

########## Data ##########
# Type            Nuclear |    Combined Cycle   |         Coal        |     Cogeneration    |       Fuel/Gas      |        Biomass
N           = [         7,         6,        46,         3,         3,       300,       300,        40,        20,        10,         1 ] #        - Number of units
C_m         = [        10,        90,        90,       100,       100,        58,        58,       150,       150,        90,        90 ] # €/MWh  - Marginal cost
C_nl        = [         0,       800,      1000,      1000,      1500,       100,       150,       300,       800,       500,      1000 ] # €/h    - No-load cost
Pg_Gen_lb   = [       950,        60,       120,       100,       350,         1,        10,        10,        30,        10,       100 ] # MW     - Power lower bound
Pg_Gen_ub   = [      1020,       300,       450,       200,       450,        10,        50,        30,        50,        25,       100 ] # MW     - Power upper bound
Pg_rr       = [       180,       600,      1200,       200,       550,        60,       120,        60,       150,        60,       250 ] # MW/h   - Maximum ramp rate for each generator of the cluster
Hg          = [         8,         5,         5,         5,         5,         4,         4,         4,         4,         3,         4 ]# s      - Inertia constant
T_st        = [      1000,         1,         2,        10,        10,         1,         1,         1,         1,         2,         4 ]# h      - Startup time
T_mut       = [        24,         2,         2,         4,         6,         1,         1,         1,         1,         2,         4 ]# h      - Minimum up time
T_mdt       = [        24,         1,         2,         2,         4,         1,         1,         1,         1,         2,         4 ]# h      - Minimum down time

Rg_max      = Pg_Gen_ub.*0.05   # MW    - FR provision
Tg          = 8                 # s     - FR delivery time

# Hora:         00:00   01:00   02:00   03:00   04:00   05:00   06:00   07:00   08:00   09:00   10:00   11:00   12:00   13:00   14:00   15:00   16:00   17:00   18:00   19:00   20:00   21:00   22:00   23:00
# Total of power demand
Pd          = [ 32.19,  28.98,  26.33,  25.04,  24.29,  22.83,  21.35,  20.40,  19.99,  20.09,  20.49,  20.38,  20.53,  21.40,  22.74,  23.66,  24.38,  25.33,  24.95,  24.25,  24.27,  26.33,  28.43,  29.77] .*10^3 # MW
# Capacity factor of RES
cf_RES      = [  6.97,   6.97,   7.08,   7.12,   6.68,   7.05,   6.81,   6.82,   6.99,   8.40,  12.94,  17.74,  19.66,  19.81,  19.85,  19.50,  18.66,  14.73,  10.06,  10.44,  10.69,   8.15,   6.82,   6.98] ./ 100

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
if (length(C_m) == length(C_nl) == length(Pg_Gen_lb) == length(Pg_Gen_ub) == length(Rg_max) == length(Hg))
    gTypes = length(C_m)
else
    println("ERROR: Wrong data input\n")
end


########## Model creation ##########
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "OutputFlag", 0)


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

        # println("RES power supply = ", round(value(P_RES * cf_RES[t] - P_curt[t]), digits = 3), " MW")
        # println("RES accommodated = ", round(value(P_curt[t]), digits = 2), " MW")
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
    plot(x, y[1, :], ylims=(0, maximum(Pd/1000)), fillrange=0, lw=0, label="Cluster 1", xlabel="Time [h]", ylabel="Power Output (MW)")
    for g in 2:gTypes
        y[g, :] .= y[g-1, :] .+ Pg_values[g, :]
        plot!(x, y[g, :], fillrange=y[g-1, :], lw=0, label="Cluster $g")
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