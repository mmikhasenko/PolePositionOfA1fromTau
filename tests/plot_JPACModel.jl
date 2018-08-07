### test
push!(LOAD_PATH, "src")
using JPACModel
using Plots
theme(:juno)

plot(x->real(t_symm(x,1.2,4)),0.2,3)
plot!(x->imag(t_symm(x,1.2,4)),0.2,3)



@time let sxv = linspace(0.3, 2.0, 60), syv = linspace(-0.25, 0.1, 60)
    cal = [qtb(sx+1im*sy, -1) for sy in syv, sx in sxv]
    contour(sxv,syv, log.(abs.(cal)), levels=50)
end

@time let sxv = linspace(0.0,1.0,40), sy = -0.25
    # cal = [interference(sx+1im*sy, 1) for sx in sxv]
    # plot(sxv, real.(cal))
    # cal = [interference(sx+1im*sy, 3) for sx in sxv]
    # plot!(sxv, real.(cal))
    cal = [interference(sx+1im*sy, 5) for sx in sxv]
    plot(sxv, real.(cal))
    # cal = [interference(sx+1im*sy, -0.1) for sx in sxv]
    # plot!(sxv, real.(cal))
end

# @time let sx = 1.1, syv = linspace(-1.0, 0.1, 20)
@time let sx = 1.3, syv = linspace(-0.8, 0.1, 20)
    cal = [interference(sx+1im*sy) for sy in syv]
    plot(syv, real.(cal), m=:c)
    # scatter!(syv, real.(cal))
    # cal = [interference(sx+1im*sy, 1.0) for sy in syv]
    # plot!(syv, real.(cal))
    # cal = [interference(sx+1im*sy, 5) for sy in syv]
    # plot!(syv, real.(cal))
end


let arr = interference_array(1.3-0.6im, 1, 150)
    N = size(arr, 1)
    Nh = Integer(N / 2)
    Nhb, Nha = Nh + [-1,+1]
    plot(
        heatmap(real.(arr[Nha:end,1:Nhb]), c=:heat, xlab="s1", ylab="s3", colorbar=false),
        heatmap(real.(arr[Nha:end,Nha:end]), c=:heat, xlab="s1", ylab="s3", colorbar=false),
        heatmap(real.(arr[1:Nhb,1:Nhb]), c=:heat, xlab="s1", ylab="s3", colorbar=false),
        heatmap(real.(arr[1:Nhb,Nha:end]), c=:heat, xlab="s1", ylab="s3", colorbar=false))
    # heatmap(imag.(arr), c=:heat)
    # real.(s1[1,:])
end


@profile interference(1.3+1im*0.0)

interference(1.3)-interference(1.3+1im*0.0)
# interference(1.3+1im*0.0)-interference(1.3-0.000001im)
interference(1.3, 1)-interference(1.3)
interference(1.3, 1)-interference(1.3-1e-5im)

plot(interference.(linspace(0.5,1.2,10)))
plot!(real.(interference.(linspace(0.5,1.2+0.0im,10))))

interference(1.1)
@time interference(1.1-0.4im)
@profile interference(1.1-0.4im)
@time interference(1.1-0.4im)

Juno.profiler()
Juno.Profile.clear()

plot!(s->real(t_symm(s+1e-5im,1.2,4)),0.2,3)
plot!(x->imag(t_symm(x+1e-5im,1.2,4)),0.2,3)

plot(x->real(t_symm_disp(x,1.2,4)),0,0.5)
plot!(x->imag(t_symm_disp(x,1.2,4)),0,0.5)
plot!(x->real(t_symm_disp(x+1e-5im,1.2,4)),0,0.5)
plot!(x->imag(t_symm_disp(x+1e-5im,1.2,4)),0,0.5)

@time let sxv = linspace(0.3, 3.0, 10), syv = linspace(-1.0, 0.1, 10)
    cal = [sy>0 ? t_symm_disp(sx+1im*sy,1.2,7) : t_symm_disp_II(sx+1im*sy,1.2,4) for sy in syv, sx in sxv]
    contour(sxv,syv, log.(abs.(cal)), levels=50)
end
