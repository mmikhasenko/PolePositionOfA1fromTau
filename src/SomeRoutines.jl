module SomeRoutines

using Juno
using Interpolations

export make_lookup_function, make_lookup_with_asymptotics_and_tan_map
export calculate_on_grid, get_results
export build_function_from_lookup
# ----------------------------------------------------------------------
# -------------------------- LOOKUP TABLES -----------------------------
"""
    function make_lookup_function(heavy_function, range)

    The method returns the function which gives quick values from the precalculated table making linear interpolation.
"""
function make_lookup_function(heavy_function, range)
    cal = [heavy_function(s) for s in range]
    lookup = interpolate((range,), cal, Gridded(Linear()))
    myf(s::Real) = (range[1] <= s <= range[end]) ? lookup[s] : error("the requested value ($s) is outside of the lookup range [$(range[1]), $(range[end])]")
    myf
end


"""
    function make_ρlookup(heavy_function, sth, lim, N=100, asympt=1.0/(8π), tol = 0.05)

    The method returns the function which gives quick values from the precalculated table making linear interpolation.
    Below the range it gives `low_asympt`-value, above the range it gives 'high_asympt'-value.
    It also checks if the provided asymptotics make sense.
"""
function make_lookup_with_asymptotics_and_tan_map(heavy_function, sth, lim, N=100,
                        low_asympt=0.0, high_asympt=1.0/(8π), tol = 0.05)
    sv = [sth+tan(t) for t in linspace(0.0,atan(lim-sth),N)]
    lf = make_lookup_function(heavy_function, sv)
    Δrel_high = abs((heavy_function(lim)-high_asympt)/high_asympt)
    Δrel_low  = abs((heavy_function(sth)-low_asympt)/ low_asympt)
    ((Δrel_high > tol) || (Δrel_low > tol)) && error("the function did not reach asymptotics, Δrels = ($(Δrel_low), $(Δrel_high)). Check the limits!")
    function myf(s::Real)
        (s < sv[1]) && return low_asympt;
        (s > sv[end]) && return high_asympt;
        return lf(s)
    end
    myf
end


# -----------------------------------------------------------------------------------
"""
    function calculate_on_grid(func, grid_x, grid_y, output)

    The method calculate 2D function on a grid displaying progress.
    The output is written to 4 text files:
      * 'output*"_re.txt"' : matrix for the real part,
      * 'output*"_im.txt"' : matrix for the imag part,
      * 'output*"_xv.txt"' : vector of the grid X-coordinates,
      * 'output*"_yv.txt"' : vector of the grid Y-coordinates.

"""
function calculate_on_grid(func, grid_x, grid_y, output)
    # interference with round 1.0
    sv = [sx+1im*sy for sy in grid_y, sx in grid_x]
    sv = hcat(sv...)
    vals = Array{Complex{Float64}}(length(sv))
    @progress for (i,s) in enumerate(sv)
        vals[i] = func(s);
    end
    cal = reshape(vals, length(grid_y), length(grid_x))
    output_re, output_im = output*"_re.txt", output*"_im.txt"
    output_xv, output_yv = output*"_xv.txt", output*"_yv.txt"
    writedlm(output_re, real.(cal))
    writedlm(output_im, imag.(cal))
    writedlm(output_xv, collect(grid_x))
    writedlm(output_yv, collect(grid_y))
    print("Calculation is complete: $(length(grid_x)*length(grid_y)) points. ")
    println("The files $output_re, $output_im, $output_xv, $output_yv are created.")
end

"""
    function get_results(output)

    The function reads result of the previous method from the file and return
    complex matrix, and the grid coordinates (X then Y)

"""
function get_results(output)
    output_re, output_im = output*"_re.txt", output*"_im.txt"
    output_xv, output_yv = output*"_xv.txt", output*"_yv.txt"
    cal = readdlm(output_re) + 1im*readdlm(output_im)
    grid_x = readdlm(output_xv)[:,1]
    grid_y = readdlm(output_yv)[:,1]
    cal, grid_x, grid_y
end


"""
    build_function_from_lookup(output)

    the function creates 2d interpolated function on the grid taking the input from the text files.

"""
function build_function_from_lookup(output)
    (cal, grid_x, grid_y) = get_results(output)
    intrp = interpolate((grid_x, grid_y), cal.', Gridded(Linear()))
    function fun(s)
        x = real(s)
        y = imag(s)
        (grid_x[1] ≤ x ≤ grid_x[end]) && (grid_y[1] ≤ y ≤ grid_x[end]) && return intrp[x,y]
        return -42.0 #error("Outside of precalculated range!")
    end
    fun
end


end
