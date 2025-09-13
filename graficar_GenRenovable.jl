# Packages
using Plots, Measures
# https://demanda.ree.es/visiona/peninsula/demandaau/acumulada/2025-01-20
# Hora:          00:00    01:00    02:00    03:00    04:00    05:00    06:00    07:00    08:00    09:00    10:00    11:00    12:00    13:00    14:00    15:00    16:00    17:00    18:00    19:00    20:00    21:00    22:00    23:00
cf_RES      = [ 11.694,  11.340,  10.160,   9.449,   8.701,   8.874,  10.256,  12.083,  13.416,  15.781,  22.351,  27.494,  28.900,  29.133,  29.503,  30.594,  30.308,  25.873,  22.033,  22.490,  23.513,  24.200,  24.049,  23.598] ./ 100
# Power from Renewable Energy Source (RES)
P_RES       = 86 # GW (Max power installed)

t = 0:23

y = cf_RES .* P_RES
# Plot
plot(t, y,
    xlabel = "Hora [h]",
    ylabel = "Potencia renovable [GW]",
    title = "Producción de RES a lo largo del día",
    legend = false,
    linewidth = 6,
    left_margin=8mm,
    bottom_margin=6mm,
    ylims=(0, 28),
    size=(1100, 650))