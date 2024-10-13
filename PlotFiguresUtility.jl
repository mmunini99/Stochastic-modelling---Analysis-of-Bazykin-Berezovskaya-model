module PlotFiguresUtility


export PlotIndividualTrajectory, PlotPredPreyTrajectories, PlotPhasePortraits, PlotEndingTrajectories, plot_simulative_pdf_3D, plot_probabilities_exiting, heatmap_link_analysis_lm

import Plots as PL
PL.plotlyjs()

function PlotIndividualTrajectory(time_line::Vector{Float64}, traj_array::Vector{Float64}, new_title::String)

    PL.plot(time_line, traj_array, label="", xlabel="Time", ylabel="Population", title=new_title)

end

function PlotPredPreyTrajectories(time_line::Vector{Float64}, traj_prey::Vector{Float64}, traj_pred::Vector{Float64}, new_title::String)

    PL.plot(time_line, traj_prey, label="Prey", legend=:topright, legend_alpha=0.01)
    PL.plot!(time_line, traj_pred, label="Predator")
    PL.xlabel!("Time")
    PL.ylabel!("Population")
    PL.title!(new_title, titlefont=PL.font(10))
    

end


function PlotPhasePortraits(tensor,  new_title::String)

    @assert length(size(tensor)) == 3 # control is 3D

    shape_tensor_3d = size(tensor)


    p = PL.plot(tensor[:,2,1], tensor[:,3,1], label="", linecolor=:blue)
    for ind in 2:shape_tensor_3d[3]
        PL.plot!(p, tensor[:,2,ind], tensor[:,3,ind], label="", linecolor=:blue)
    end
    PL.xlabel!("Prey")
    PL.ylabel!("Predator")
    PL.title!(new_title, titlefont=PL.font(10))
    PL.display(p)
    
end



function PlotEndingTrajectories(tensor_set, x_vector::Vector{Float64}, dim::String, new_title::String, name_x::String, name_y::String, col_vec)

    @assert all(ndims(tensor) == 3 for tensor in tensor_set) # control is 3D

    @assert dim in ["prey", "predator"] "Choice is prey or predator"

    

    @assert all(size(tensor)[3] == length(x_vector) for tensor in tensor_set)

    dim_ind = dim == "prey" ? 1 : 2

    p = PL.scatter(size=(1400, 600), markersize=3)


    for tensor_ind in 1:length(tensor_set)


        dim_end_obs = size(tensor_set[tensor_ind])[1]
        x_data = repeat(x_vector, inner=dim_end_obs)
        y_values = vec(tensor_set[tensor_ind][:,dim_ind,:])
        PL.scatter!(p, x_data, y_values, label="", markercolor=col_vec[tensor_ind], alpha = 0.1)



    end

    PL.xlabel!(name_x)
    PL.ylabel!(name_y)
    PL.title!(new_title, titlefont=PL.font(10))


    PL.display(p)


    
end


function plot_simulative_pdf_3D(tensor, title::String)

    first = sort(unique(vec(tensor)), rev=true)[1]
    second = sort(unique(vec(tensor)), rev=true)[2]
    new_first = first/second > 4 ? second * 2 : first
    tensor[1,1] = new_first

    x = collect(range(0, stop=1, length=size(tensor)[1]))
    pl = PL.heatmap(x, x, tensor, c=PL.cgrad(:heat),  size=(600, 600), colorbar = false)
    PL.ylabel!("prey")
    PL.xlabel!("predator")
    PL.title!("heat", titlefont=PL.font(10))
    PL.display(pl)
    
end

function plot_probabilities_exiting(matrix::Array{Float64, 2}, vector::Vector{Float64}, title::String)
    pl = PL.heatmap(vector, vector, matrix, c=PL.cgrad(:heat),  size=(600, 600), colorbar = false)
    for i in 1:length(vector)
        for j in 1:length(vector)
            PL.annotate!(pl, ((i-1)/10 + 0.05, (j-1)/10 + 0.05), PL.text(string(round(matrix[i, j], digits=2)), 8, :black))
        end
    end
    PL.xlabel!("Magnitude of l's noise")
    PL.ylabel!("Magnitude of m's noise")
    PL.title!(title, titlefont=PL.font(10))
    PL.display(pl)
end


function heatmap_link_analysis_lm(matrix, l_vector, m_vector)

    sets_names = String[] 

    for ind_l in l_vector
        for ind_m in m_vector
            l = round(ind_l, digits=2)
            m = round(ind_m, digits=2)
            push!(sets_names, "$l-$m")
        end
    end

    PL.heatmap(sets_names, sets_names, matrix, size=(1400, 1000))


end



end