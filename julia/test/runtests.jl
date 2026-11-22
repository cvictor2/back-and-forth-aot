using Test, BFNPaper, MAT, JSON3, LinearAlgebra
const B=BFNPaper
@testset "MATLAB reference trajectories" begin
    cases=B.as_records(matread(joinpath(@__DIR__,"fixtures","matlab_short_runs.mat"))["cases"])
    @test length(cases)==31
    for f in cases
        @testset "$(f["name"])" begin
            c=B.config_from_dict(JSON3.read(f["config_json"]));c.save_output=false;c.make_plots=false
            r=run_experiment(c)
            @test r["status"]=="completed"
            for (a,b) in (("initial_reference","initial"),("terminal_reference","terminal"),("final_state","final"))
                @test vec(r[a]) ≈ vec(f[b]) rtol=2e-11 atol=2e-10
            end
            @test r["diagnostics"][:,6:11] ≈ f["errors"] rtol=2e-10 atol=2e-11
            @test r["cycle_errors"] ≈ f["endpoints"] rtol=2e-10 atol=2e-11
        end
    end
end
@testset "Independent numerical properties" begin
    c=preset(:heat,:high);c.n=128;c.nonlinear_cutoff=42;c.high_mode=30;c.viscosity_cutoff=20
    c.dt=1e-5;c.window=.0002;c.spinup=0;c.cycles=1;c.mu_forward=c.mu_backward=0;c.initialization=:exact;c.save_output=false
    g=B.spectral_grid(c);k=B.kernel(g,c,-c.dt;backward=true)
    @test k.E[31]≈1 # Only the viscosity is truncated, not this state coefficient.
    @test abs(B.spectral_step(g.u0,g.u0,g,c,k)[31])>1
    e=complex.(collect(1:c.n),collect(-c.n:-1));d=B.error_values(e,g,c)
    @test d[1]^2≈d[2]^2+d[3]^2
    @test d[4]^2≈d[5]^2+d[6]^2
    @test d[1]≈sqrt(c.length)*norm(e)/c.n
    c.scheme=:voigt;c.alpha=.1;k=B.kernel(g,c,-c.dt;backward=true)
    @test k.V≈1 ./ (1 .+ c.alpha^2 .*g.k2)
    c.filtered=true;k=B.kernel(g,c,-c.dt;backward=true)
    @test all(k.V[g.obs].==1)
    c=preset(:lorenz);c.schedule=:thirds;c.components=[1,2,3];c.fractions=[0,0,0]
    @test !any(B.observation_mask(c,0.0))
    c.fractions=[1,1,1]
    @test B.observation_mask(c,0.0)==[true,false,false]
    @test B.observation_mask(c,1/3)==[false,true,false]
    @test !any(B.observation_mask(c,1.0))
    c=preset(:nse);c.n=64;c.observation_cutoff=12;c.nonlinear_cutoff=21;c.allow_restart_resample=true
    g=B.spectral_grid(c);x=(0:63).*2pi/64;z=g.forward_plan*complex.(cos.(x).*transpose(cos.(x)))
    d=B.error_values(z,g,c)
    @test d[1]≈sqrt(2)*pi
    @test d[4]≈2pi
    @test !any(g.dx[33,:].!=0) && !any(g.dy[:,33].!=0)
    # Nonzero analytic advection independently checks the streamfunction sign.
    psi=g.forward_plan*complex.(cos.(x) .+ 0.2 .*transpose(cos.(2 .*x)))
    q=B.nonlinear(psi,g,c)-g.forcing
    expected=-0.24 .*sin.(x).*transpose(sin.(2 .*x))
    @test real.(g.inverse_plan*q) ≈ expected atol=1e-12
    u=-real.(g.inverse_plan*(g.dy.*psi));v=real.(g.inverse_plan*(g.dx.*psi))
    omega=-g.k2.*psi
    advection=u.*real.(g.inverse_plan*(g.dx.*omega))+v.*real.(g.inverse_plan*(g.dy.*omega))
    @test q ≈ -g.invlap.*(g.forward_plan*advection) rtol=1e-12 atol=1e-11

    @test_throws ErrorException B.config_from_dict(Dict("not_an_option"=>1))
    c=preset(:transport);c.dt=0;@test_throws ArgumentError B.validate(c)
    c=preset(:transport);c.max_reference_bytes=1;c.save_output=false;@test_throws ErrorException run_experiment(c)
    c=preset(:kdv_damped);c.n=128;c.nonlinear_cutoff=42;c.window=.0002;c.cycles=1;c.max_amplitude=.01;c.save_output=false
    r=run_experiment(c);@test r["status"]=="failed";@test size(r["cycle_errors"],1)==0
end
@testset "Portable result round trip" begin
    mktempdir() do dir
        c=preset(:transport);c.n=128;c.nonlinear_cutoff=42;c.window=.0002;c.dt=1e-5;c.spinup=0;c.cycles=1;c.make_plots=false;c.output_root=dir
        r=run_experiment(c);q=matread(joinpath(r["output_folder"],"result.mat"))["result"]
        @test q["final_state"]==r["final_state"]
        @test q["diagnostics"]==r["diagnostics"]
        @test B.config_from_dict(JSON3.read(q["config_json"])).model==:transport
        @test isfile(joinpath(r["output_folder"],"cycle_errors.csv"))
    end
end

@testset "Output survives process exit" begin
    # Inspect output only after the producing process exits. This verifies
    # persistence separately from the in-process MAT round-trip checks.
    mktempdir() do parent
        script=joinpath(parent,"producer.jl");marker=joinpath(parent,"folder.txt")
        write(script,"""
        using BFNPaper
        suite=BFNPaper.output_directory($(repr(parent)),"suite_")
        c=preset(:transport);c.n=128;c.nonlinear_cutoff=42;c.dt=1e-5
        c.window=.0002;c.spinup=0;c.cycles=1;c.make_plots=false;c.output_root=suite
        r=run_experiment(c)
        write($(repr(marker)),r["output_folder"])
        """)
        project=dirname(@__DIR__)
        run(`$(Base.julia_cmd()) --startup-file=no --project=$project $script`)
        folder=read(marker,String)
        @test isdir(dirname(folder))
        @test isfile(joinpath(folder,"result.mat"))
        @test isfile(joinpath(folder,"cycle_errors.csv"))
        saved=matread(joinpath(folder,"result.mat"))["result"]
        @test saved["completed_cycles"]==1
        @test all(isfinite,saved["final_state"])
    end
end
