push!(LOAD_PATH, "src")
using JPACModel
using Plots
theme(:juno)

@time interference(1.3-0.3im,  1.0)
@time interference(1.3-0.3im,  3.0)

using SomeRoutines

calculate_on_grid(s->interference(s,5),
                  linspace(0.3, 2.0, 2), linspace(-1.0, 0.1, 3),
                  "test")

# interference with round 1.0
calculate_on_grid(s->interference(s,1),
                  linspace(0.3, 2.0, 100), linspace(-1.0, 0.1, 100),
                  "new_interf_rr1")

# interference with round 3.0
calculate_on_grid(s->interference(s,3),
                  linspace(0.3, 2.0, 100), linspace(-1.0, 0.1, 100),
                  "new_interf_rr3")

# interference with round 5.0
calculate_on_grid(s->interference(s,5),
                  linspace(0.3, 2.0, 100), linspace(-1.0, 0.1, 100),
                  "new_interf_rr5")

# qtb with round 3.0
calculate_on_grid(s->qtb(s,3),
                  linspace(0.3, 2.0, 100), linspace(-1.0, 0.1, 100),
                  "qtb_rr3")

# qtb with round 3.0
calculate_on_grid(s->qtb(s,5),
            linspace(0.3, 2.0, 100), linspace(-1.2, 0.1, 130),
            "qtb_rr5_rect")

let (cal, grid_x, grid_y) = get_results("interf_rr3")
    contour(grid_x, grid_y, log.(abs.(cal)), levels=50)
end

let (cal, grid_x, grid_y) = get_results("qtb_rr5_rect")
    contour(grid_x, grid_y, log.(abs.(cal)), levels=50)
    # grid_x, grid_y
end

let (cal, grid_x, grid_y) = get_results("new_interf_rr5")
    contour(grid_x, grid_y, log.(abs.(cal)), levels=50)
end

let (cal, grid_x, grid_y) = get_results("new_interf_rr5")
    contour(grid_x, grid_y, log.(abs.(cal)), levels=50)
    yInd = 60;
    y = grid_y[yInd];
    plot( grid_x, log.(abs.(cal[yInd,:])))
    plot!(grid_x, log.(abs.([interference(x+1im*y, 5) for x in grid_x])))
end

let (cal, grid_x, grid_y) = get_results("interf_rr3")
    contour(grid_x, grid_y, log.(abs.(cal)), levels=50)
end

let (cal, grid_x, grid_y) = get_results("new_interf_rr1")
    N = 90
    cal1 = cal[N,:]
    plot(grid_x, real.(cal1), lab="si = $(grid_y[N])")
    plot!(grid_x, imag.(cal1), lab="si = $(grid_y[N])")
end

let (cal, grid_x, grid_y) = get_results("interf_rr1")
    N = 90
    cal1 = cal[N,:]
    plot!(grid_x, real.(cal1), lab="si = $(grid_y[N])")
    plot!(grid_x, imag.(cal1), lab="si = $(grid_y[N])")
end


let (cal, grid_x, grid_y) = get_results("qtb_rr5")
    contour(grid_x, grid_y, log.(abs.(cal)), levels=50)
end

qtbi = build_function_from_lookup("qtb_rr3")
plot(x->real(qtbi(x-0.25im)), 0.3, 2.0)
plot!(x->imag(qtbi(x-0.25im)), 0.3, 2.0)

interfi = build_function_from_lookup("interf_rr3")
plot!(x->real(interfi(x-0.25im)), 0.3, 2.0)
plot!(x->imag(interfi(x-0.25im)), 0.3, 2.0)
