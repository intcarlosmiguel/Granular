using Pkg
using DataFrames, CSV, DelimitedFiles
using Plots

# Carregar os dados do arquivo example.txt
data = readdlm("./example.txt")'
print("Acabou!")
df = DataFrame(Tempo=data[1, :], id=Int.(data[2, :]), x=data[3, :], y=data[4, :])
# Configurações iniciais do plot
dt = 0.001
cores = ["red", "blue", "green"]
m = 150 / 2
angulo = 45
p = plot(size=(800, 800), xlim=(-324 + m, 324 + m), ylim=(0, 910), aspect_ratio=:equal)

# Elementos estáticos
plot!([0, 0], [0, 98], color=:blue)
plot!([150, 150], [0, 98], color=:blue)
plot!([0, -154], [98, 98 + 154 * tand(angulo)], color=:blue)
plot!([150, 304], [98, 98 + 154 * tand(angulo)], color=:blue)
plot!([-154, -154], [98 + 154 * tand(angulo), 910], color=:blue)
plot!([304, 304], [98 + 154 * tand(angulo), 910], color=:blue)

# Preparar as partículas
particulas = [scatter!([], [], marker=(:circle, 7.5), color=:blue) for _ in 1:nrow(df)]

# Função de atualização para animação
function update(t)
    filtro = df[df.Tempo .== t, :]
    for i in eachrow(filtro)
        particulas[i.id][1][1] = i.x * 1000
        particulas[i.id][1][2] = i.y * 1000
        set_marker_color(particulas[i.id], cores[(i.id % 3) + 1])
    end
end


# Criar animação
anim = @animate for t in minimum(df.Tempo):dt:maximum(df.Tempo)
    update(t)
    frame(p)
end

# Salvar como GIF
gif(anim, "animation.gif", fps=200)
