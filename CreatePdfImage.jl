include("C:\\Users\\matte\\Documents\\STO MOD AND SIM\\Progetto Esame\\SolutionSystem.jl")
import .SimulationStochasticEnding as SSEnd

using Makie
using GLMakie
using Images



function CreatePdfImage(r, l, k, m, param1, param2)
    Z = SSEnd.normalized_simulative_pdf(r, l, k, m, param1, param2, 1000.0, 10000, [0.0, 1.0], 10000, 100)

    println(pwd()) 
    

    first = sort(unique(vec(Z)), rev=true)[1]
    second = sort(unique(vec(Z)), rev=true)[2]
    if second > 0.0
        new_first = first/second > 4 ? second * 2 : first
        Z[1,1] = new_first
    end

    fig = Makie.Figure(size = (600, 600))

    # Create an axis within the figure
    ax = Makie.Axis(fig[1, 1], title = "Noise l: $param1, Noise m: $param2", xlabel="prey", ylabel="predator")


    # Plot the heatmap on the axis
    Makie.heatmap!(ax, Z, colormap = :viridis)

    png_file = "frame_r$(round(r, digits=2))_l$(round(l, digits=2))_k$(round(k, digits=2))_m$(round(m, digits=2))_frame_$(round(param1, digits=2))_$(round(param2, digits=2)).png"
    save(png_file, fig)

    return 0

end