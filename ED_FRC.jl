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
Rg_max      = [   0,  200,    0] # MW       - FR provision
Hg          = [   5,    5,    5] # s        - Inertia constant
Tg          = 10                 # s        - FR delivery time

# Frecuency data
f_0         = 50 # Hz
Δf_max      = 0.5 # Hz
Δfss_max    = 0.2 # Hz

# Power from Renewable Energy Source (RES)
P_RES       = 40 *10^3 # MW (Max power installed)
cf_RES      = 0.25 # Capacity factor of RES (<= 1)
Rs_ub       = 0 # MW (EFR)
Ts          = 1 # s

# Maximum power infeed
P_L_max = maximum(Pg_Gen_ub)

# Load damping
D = 1.5 # %/Hz

# Total of power demand
Pd = 24 *10^3 # MW


# Number of total generators
G   = sum(N)


# Check if the length of all data are the same
if (length(Pg_cost) == length(Pg_nl_cost) == length(Pg_Gen_lb) == length(Pg_Gen_ub) == length(Rg_max) == length(Hg))
    aux = length(Pg_cost)
else
    println("ERROR: Wrong data input")
end

# Generation of arrays
Pg_cost = vcat([fill(Pg_cost[i], N[i]) for i in 1:aux]...)
Pg_nl_cost = vcat([fill(Pg_nl_cost[i], N[i]) for i in 1:aux]...)
Pg_Gen_lb = vcat([fill(Pg_Gen_lb[i], N[i]) for i in 1:aux]...)
Pg_Gen_ub = vcat([fill(Pg_Gen_ub[i], N[i]) for i in 1:aux]...)
Rg_max = vcat([fill(Rg_max[i], N[i]) for i in 1:aux]...)
Hg = vcat([fill(Hg[i], N[i]) for i in 1:aux]...)


########## Model creation ##########
model = Model(Gurobi.Optimizer)
set_silent(model)


########## Variables ##########
@variable(model, y[1:G], Bin) # on/off of each generator
@variable(model, Pg[1:G]) # Power from each generator
@variable(model, Rg[1:G]) # PRF from each generators
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
@objective(model, Min, sum(y[n] * Pg_nl_cost[n] + Pg[n] * Pg_cost[n] for n in 1:G))


########## Constraints ##########
# @constraint(model, sum(Pg) + P_RES - P_curt == Pd - loadShedding)
@constraint(model, sum(Pg) + P_RES * cf_RES - P_curt == Pd)

# # Number of active generators of each group constraint
# @constraint(model, [i in 1:G], 0 <= Ng[i])
# @constraint(model, [i in 1:G], Ng[i] <= Ng[i])

# Power from each generator constraints
@constraint(model, [i in 1:G], Pg_Gen_lb[i] * y[i] <= Pg[i])
@constraint(model, [i in 1:G], Pg[i] <= Pg_Gen_ub[i] * y[i])

# PRF provision from g constraints
@constraint(model, [i in 1:G], 0 <= Rg[i])
@constraint(model, [i in 1:G], Rg[i] <= y[i] * Rg_max[i])
@constraint(model, [i in 1:G], Rg[i] <= Pg_Gen_ub[i] - Pg[i])

# EFR provision from RES constraints
@constraint(model, 0 <= P_curt)
@constraint(model, P_curt <= P_RES * cf_RES)
@constraint(model, 0 <= Rs) # Rs (storage)
@constraint(model, Rs <= Rs_ub)


# Largest power infeed constraints
@constraint(model, P_L <= P_L_max)
for i in 1:G
    @constraint(model, P_L >= Pg[i])
end


# Equation 8: System inertia
@constraint(model, H == sum(Hg .* Pg_Gen_ub .* y))

# Equation 9
# @constraint(model, (P_L - Rs - sum(Rg) / (D * Pd)) <= Δfss_max)
@constraint(model, Rs + sum(Rg) >= P_L - Δfss_max * D * Pd)

# Equiation 13
@constraint(model, (H / f_0 - Rs * Ts / (4*Δf_max)) * sum(Rg) >= (P_L - Rs)^2 * Tg / (4*Δf_max))

# Equation 19
# @constraint(model, (H / f_0 - Rs * Ts / (4*Δf_max)) * sum(Rg) >= (P_L - Rs)^2 * Tg / (4*Δf_max) - (P_L - Rs) * Tg * D * Pd / 4)


########## Optimization ##########
# relax_integrality(model)
optimize!(model)


########## Solution display ##########
clearTerminal()

if termination_status(model) == OPTIMAL || termination_status(model) == LOCALLY_SOLVED
    println("##### Solución encontrada #####")
    println(termination_status(model))
    println("H / f_0 - Rs * Ts / (4*Δf_max) = ", round(value((Hg' * Pg_Gen_ub) / f_0 - Rs * Ts / (4*Δf_max)), digits = 3))

    # for i in 1:G
    #     println("Rg del grupo de generador ", i, " = ", round(value(Rg[i]),digits = 3))
    # end
    println("RG = ", round(value(sum(Rg)), digits = 3))

    println("\nDemanda = ", Pd)

    println("Potencia total generada = ", round(sum(value(Pg[i]) for i in 1:G), digits = 3))

    println("Potencia RES = ", round(value(P_RES * cf_RES - P_curt), digits = 3))
    println("P_curt = ", round(value(P_curt), digits = 3))

    for i in 1:G
        println("Potencia generador ", i, " : ", round(value(Pg[i]), digits = 3), " MW")
    end

    print(Pg_cost)

    idx = cumsum([1; Ng[1:end-1]])
    for i in 1:length(N)
        gen_idx = idx[i]:idx[i] + Ng[i] - 1
        println("Nomber of generators of type ", i, " = ", sum(value(y[n]) for n in gen_idx))
    end

    println("Total cost: ", round(objective_value(model), digits = 3))

# elseif value(H / f_0 - Rs * Ts / (4*Δf_max)) < 0
    # println("ERROR: La combinación de inercia 'H' y la respuesta rápida 'Rs' no es suficiente para cubrir la desviación de frecuencia")

else
    println("ERROR: ", termination_status(model))

end