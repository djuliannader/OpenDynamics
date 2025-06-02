module plotwigner

using Plots
export plotingw, plot_function

# Function to read data from a file
function readdata(filename::String)
    x = Float64[]
    p = Float64[]
    Q = Float64[]
    open(filename) do f
        for line in eachline(f)
            ssplited = split(line)
            push!(x, parse(Float64, ssplited[1]))
            push!(p, parse(Float64, ssplited[2]))
            push!(Q, parse(Float64, ssplited[3]))
        end
    end
    return [x, p, Q]
end

# Function to plot the combined heatmap and contour
function plot_combined(input_name::String, Delta::Float64, k::Float64, e1::Float64, output_name::String, title::String)
    # Read data for the heatmap
    a = readdata(input_name)
    nn = floor(Int, sqrt(length(a[1])))
    xx = reshape(a[1], (nn, nn))
    pp = reshape(a[2], (nn, nn))
    ww = reshape(a[3], (nn, nn))
    x_vals = unique(a[1])
    p_vals = unique(a[2])

    # Define the contour function
    f(x, y) = -Delta/2 * (x^2 + y^2) + k/4 * (x^2 + y^2)^2 - e1 * (x^2 - y^2) - 24

    # Create the heatmap and overlay the contour plot
  #  heatmap(x_vals, p_vals, ww, color=:bwr, clims=(-0.1, 0.1), 
  #          aspect_ratio=1, title=title, xlims=(-8, 8), ylims=(-8, 8), 
  #          framestyle=:box)          with title
    heatmap(x_vals, p_vals, ww, color=:bwr, clims=(-0.1, 0.1), 
            aspect_ratio=1,  xlims=(-8, 8), ylims=(-8, 8), 
            framestyle=:box)
   # contour!(x_vals, p_vals, f, levels=[0], linewidth=2, color=:green) Activarla para graficar la trayectoria clasica

    # Save the combined plot
    savefig(output_name)
    println("Combined heatmap and contour plot saved to ", output_name)
end

end # module

# Example usage
using .plotwigner

for i in 3.0:0.1:3.3
    # Input file and output file name
    input_file = string("output/gsDelta=", i, ".dat")
    output_combined = string("output/D = ", i, ".png")

    # Title for the combined plot
    #t_value = round(i * 0.1, digits=2)
    #title = "t = $t_value"

    Delta = -2.0
    k = 1.0
    e1 =0.0

    plotwigner.plot_combined(input_file, Delta, k, e1, output_combined, "i")
end
