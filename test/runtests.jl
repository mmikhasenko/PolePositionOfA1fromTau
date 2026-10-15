using Test, QuadGK, LinearAlgebra, DelimitedFiles
using PolePositionOfA1fromTau
using PolePositionOfA1fromTau: Isobars, Masses, JPACModel, SomeRoutines, Fitting, DalitzPlotAnalysis, ErrorPropagation

@testset "Isobar normalization and real phase space" begin
    @test qtb(0.0) == qtb(9Masses.mπ2) == 0.0
    @test interference(0.0) == 0.0
    @test quadgk(u->real(Isobars.U1(u))/(2π),4Masses.mπ2,Inf)[1] ≈ 1 rtol=1e-6
    for u in (0.2,0.6,1.1)
        @test Isobars.f1_II(u) ≈ conj(Isobars.f1_I(u))
        @test Isobars.f1_I(complex(u)) ≈ Isobars.f1_I(u)
    end
    for s in (0.3,1.0,2.0)
        @test qtb(s) > 0
        @test qtb(complex(s),0.0) ≈ qtb(s) rtol=1e-5
        @test symm(s) >= 0
        @test symm(s) ≈ qtb(s)+interference(s)
    end
    @test isfinite(interference(1.5-0.05im))
    @test size(JPACModel.interference_array(1.5-0.05im,1.0,4)) == (4,4)
    # Equal-mass Dalitz transform and angular conventions.
    s,u,z = 1.5,0.4,0.3
    v=DalitzPlotAnalysis.cross_basis_s3(u,z,Masses.mπ2,Masses.mπ2,Masses.mπ2,s)
    @test DalitzPlotAnalysis.cross_basis_z23(u,v,Masses.mπ2,Masses.mπ2,Masses.mπ2,s) ≈ z
    @test DalitzPlotAnalysis.Wignerd(1,0,0,0.7) ≈ cos(0.7)
end

@testset "Dispersive integral against an analytic logarithm" begin
    # Constant density: Pi0(s)=-log(1-s/a)/pi. Tests PV and complex integral independently.
    a=0.2; b=0.6
    model=DispersiveLineshape(x->1.0;threshold=a,subtraction=b,rtol=1e-9)
    offset=-log(abs(1-b/a))/π
    exact(s)=-log(1-complex(s)/a)/π-offset
    @test real(selfenergy(model,b)) ≈ 0 atol=1e-12
    for s in (-0.3,0.1,0.5,1.2)
        value=selfenergy(model,s)
        @test real(value) ≈ -log(abs(1-s/a))/π-offset atol=1e-9
        @test imag(value) ≈ (s>a ? 1 : 0)
    end
    for s in (0.5+0.2im,0.5-0.2im,-0.3+0.1im)
        @test selfenergy(model,s) ≈ exact(s) rtol=1e-8
        @test selfenergy(model,conj(s)) ≈ conj(selfenergy(model,s))
    end
    @test selfenergy(model,0.7+1e-5im)-selfenergy(model,0.7-1e-5im) ≈ 2im atol=2e-5
    @test selfenergy(model,0.7-0.1im;sheet=:II,continued_density=x->1) ≈ exact(0.7-0.1im)+2im
    @test_throws ArgumentError selfenergy(model,0.7-0.1im;sheet=:II)
    @test_throws ArgumentError selfenergy(model,0.7;sheet=:II,continued_density=x->1)
    @test imag(inv(amplitude(model,0.7,1.3,7))) ≈ -0.5
    @test intensity(model,0.1,1.3,7) == 0
    @test intensity(model,4.0,1.3,7) == 0
    # Continuous threshold with known integral at s=a and a genuinely varying density.
    smooth=DispersiveLineshape(x->1-a/x;threshold=a,subtraction=0.0)
    @test selfenergy(smooth,a) ≈ 1/π atol=1e-7
    @test imag(selfenergy(smooth,0.6)) ≈ 1-a/0.6
    @test real(selfenergy(smooth,0.6+1e-5im)) ≈ real(selfenergy(smooth,0.6)) atol=1e-4
end

@testset "Lookup bounds, fitting, covariance" begin
    mktempdir() do dir
        out=joinpath(dir,"grid")
        xs=[0.0,1.0,2.0];ys=[-0.5,0.0]
        SomeRoutines.calculate_on_grid(identity,xs,ys,out)
        f=SomeRoutines.build_function_from_lookup(out)
        @test f(0.3-0.2im) ≈ 0.3-0.2im
        @test_throws DomainError f(0.3+0.2im)
    end
    data=[0.0 2.0 0.1;1.0 5.0 0.1;2.0 8.0 0.1]
    chi2=Fitting.build_chi2((x,a,b)->a+b*x,data,Matrix(Diagonal(fill(0.01,3))))
    result=Fitting.fit_with_grad(chi2;starting_pars=[1.0,1.0],verbose=false)
    @test result[2] ≈ [2,3] atol=1e-4
    @test result[1] < 1e-6
    @test length(ErrorPropagation.get_ellipce_from_covariance([1.0 0;0 2.0],[0,0])) == 1000
end

@testset "Physical model boundary and legacy API agreement" begin
    model=JPACModel.dispersive_model(:qtb)
    for s in (0.5,1.2,2.0)
        @test imag(selfenergy(model,s)) ≈ model.rho(s)
        @test amplitude(model,s,1.3,7) ≈ JPACModel.t_qtb_disp(s,1.3,7)
        @test imag(inv(amplitude(model,s,1.3,7))) ≈ -model.rho(s)/2
    end
    @test real(selfenergy(model,(Masses.mπ+Masses.mρ)^2)) ≈ 0 atol=1e-9
end
