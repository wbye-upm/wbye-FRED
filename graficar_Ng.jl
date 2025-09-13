using Plots

# Eje X: horas 1–24
x = 0:24

##### Datos #####
# Carbón 1
Caso1 = [2,1,0,0,0,1,2,2,2,2,2,0,0,0,0,0,0,1,2,2,2,2,2,2,2]
Caso2 = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
Caso3 = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
Caso4 = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
N_max = 2

##### CC 2
# Caso1 = [0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
# Caso2 = [0,0,0,0,0,0,11,22,22,22,22,11,0,0,0,0,0,2,2,2,2,2,0,0,0]
# Caso3 = [5,5,5,5,5,5,9,17,22,22,17,13,5,3,1,0,0,3,6,6,6,6,3,0,0]
# Caso4 = [5,5,5,5,5,5,10,19,24,24,19,14,5,5,0,0,0,4,7,7,7,7,3,0,0]
# N_max = 33

##### CC 3
# Caso1 = [5,0,0,0,0,0,6,15,18,18,6,0,0,0,0,0,0,0,9,13,14,14,7,0,0]
# Caso2 = [5,5,5,5,4,4,8,12,14,14,10,6,2,0,0,0,0,8,17,17,17,17,16,15,15]
# Caso3 = [14,14,14,14,14,14,14,15,15,15,12,9,9,7,7,5,8,12,16,18,18,18,18,18,18]
# Caso4 = [14,14,14,14,14,14,14,14,14,14,10,9,9,5,5,1,5,9,14,18,18,18,18,18,18]
# N_max = 18

##### CC 4
# Caso1 = [20,20,19,17,18,20,20,20,20,20,20,20,15,13,13,10,13,20,20,20,20,20,20,20,20]
# Caso2 = [20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20]
# Caso3 = [20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20]
# Caso4 = [20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20]
# N_max = 20

# Subplots en cuadrícula 2x2
plt = plot(
    plot(x, Caso1, seriestype=:steppost, lw=3, label="", title="Caso 1",
         titlefont=font(11), top_margin=-1mm,
         xlabel="Tiempo [h]", ylabel="Ng",
         xticks=0:3:24, yticks=[0:5:N_max; N_max], ylims=(0,N_max+0.2), grid=true),
    plot(x, Caso2, seriestype=:steppost, lw=3, label="", title="Caso 2",
         titlefont=font(11), top_margin=-1mm,
         xlabel="Tiempo [h]", ylabel="Ng",
         xticks=0:3:24, yticks=[0:5:N_max; N_max], ylims=(0,N_max+0.2), grid=true),
    plot(x, Caso3, seriestype=:steppost, lw=3, label="", title="Caso 3",
         titlefont=font(11), top_margin=-1mm,
         xlabel="Tiempo [h]", ylabel="Ng",
         xticks=0:3:24, yticks=[0:5:N_max; N_max], ylims=(0,N_max+0.2), grid=true),
    plot(x, Caso4, seriestype=:steppost, lw=3, label="", title="Caso 4",
         titlefont=font(11), top_margin=-1mm,
         xlabel="Tiempo [h]", ylabel="Ng",
         xticks=0:3:24, yticks=[0:5:N_max; N_max], ylims=(0,N_max+0.2), grid=true),
    layout=(2,2), size=(1000,700), margin=8mm
)

display(plt)