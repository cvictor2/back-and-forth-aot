"""Explicit numerical settings. `preset` supplies the manuscript defaults."""
Base.@kwdef mutable struct Config
    name::String = "transport_mixed"
    model::Symbol = :transport
    profile::Symbol = :mixed
    n::Int = 4096
    length::Float64 = 2.0
    dt::Float64 = 1e-4
    window::Float64 = 1.0
    cycles::Int = 5
    nu::Float64 = 0.0
    gamma::Float64 = 0.0
    transport_speed::Float64 = 1.0
    dispersion::Float64 = 1.0
    integrator::Symbol = :ifeuler
    feedback_factor::Bool = false
    nonlinear_cutoff::Int = 682
    mu_forward::Float64 = 1000.0
    mu_backward::Float64 = 1000.0
    observation_cutoff::Int = 16
    full_observations::Bool = false
    components::Vector{Int} = [1, 2]
    fractions::Vector{Float64} = [1.0, 1.0, 1.0]
    schedule::Symbol = :continuous
    window_rule::Symbol = :halfopen
    backward_observation_shift::Int = 0
    high_mode::Int = 30
    high_amplitude::Float64 = 0.05
    low_amplitude::Float64 = 1.0
    spinup::Float64 = 0.0101
    initialization::Symbol = :zero
    reference_initial::Vector{Float64} = [20.0, 30.0, 40.0]
    guess_initial::Vector{Float64} = [30.0, 40.0, 50.0]
    sigma::Float64 = 10.0
    rho::Float64 = 28.0
    beta::Float64 = 8 / 3
    scheme::Symbol = :standard
    filtered::Bool = false
    alpha::Float64 = 1e-3
    viscosity_factor::Float64 = -1.0
    damping_factor::Float64 = -1.0
    viscosity_cutoff::Int = 50
    diffusion_correction::Float64 = 0.0
    cancel_damping::Bool = false
    forcing_amplitude::Float64 = 1.0
    restart_file::String = joinpath(REPO_ROOT, "data", "NSE_reference_and_forcing.mat")
    allow_restart_resample::Bool = false
    stride::Int = 10
    max_norm::Float64 = 1e100
    max_amplitude::Float64 = 1e100
    stop_tolerance::Union{Nothing,Float64} = nothing
    max_reference_bytes::Int = 2_000_000_000
    output_root::String = joinpath(REPO_ROOT, "results", "julia")
    save_output::Bool = true
    make_plots::Bool = true
    save_cycle_states::Bool = true
end

function preset(model::Symbol, profile::Symbol=:mixed)
    c = Config(model=model, profile=profile, name="$(model)_$(profile)")
    if model == :lorenz
        c.n=3; c.dt=1e-5; c.window=1; c.cycles=100; c.spinup=0
        c.integrator=:rk4; c.initialization=:explicit
        c.mu_forward=100; c.mu_backward=100; c.stride=500
    elseif model == :heat
        c.transport_speed=0; c.nu=0.05; c.scheme=:cutoff
    elseif model == :transport
        # Numerical transport direction is +a u_x, matching the MATLAB driver.
    elseif model == :burgers
        c.transport_speed=0; c.dt=1e-6; c.window=0.01; c.spinup=0.010001
    elseif model == :nse
        c.n=1024; c.length=2pi; c.dt=1e-3; c.window=0.05; c.nu=1e-4
        c.transport_speed=0; c.spinup=0; c.initialization=:projection
        c.observation_cutoff=50; c.nonlinear_cutoff=floor(Int, c.n/3); c.stride=1
    elseif model in (:kdv_damped, :kdv_viscous)
        c.n=512; c.dt=1e-5; c.window=1; c.spinup=0; c.integrator=:ifrk4
        c.feedback_factor=true; c.transport_speed=0; c.high_mode=20
        c.high_amplitude=0.1; c.observation_cutoff=15
        c.nonlinear_cutoff=floor(Int,c.n/3); c.initialization=:exact; c.stride=100
        if model == :kdv_damped
            c.gamma=0.5; c.alpha=1e-7; c.diffusion_correction=1e-5
        else
            c.nu=1e-3; c.alpha=1e-2; c.diffusion_correction=c.nu
        end
    else
        throw(ArgumentError("Unknown model $model"))
    end
    c
end

function validate(c::Config)
    c.model in (:lorenz,:heat,:transport,:burgers,:nse,:kdv_damped,:kdv_viscous) || error("Unknown model")
    for x in (c.dt,c.window,c.length,c.max_norm,c.max_amplitude)
        isfinite(x) && x>0 || throw(ArgumentError("Time, domain and stop limits must be positive and finite"))
    end
    for x in (c.nu,c.gamma,c.mu_forward,c.mu_backward,c.alpha,c.diffusion_correction)
        isfinite(x) && x>=0 || throw(ArgumentError("Coefficients must be finite and nonnegative"))
    end
    all(isfinite,(c.transport_speed,c.dispersion,c.forcing_amplitude,c.low_amplitude,c.high_amplitude,c.viscosity_factor,c.damping_factor)) || error("Nonfinite numerical coefficient")
    c.max_reference_bytes>0 || error("max_reference_bytes must be positive")
    c.stop_tolerance===nothing || (isfinite(c.stop_tolerance) && c.stop_tolerance>=0) || error("Invalid stopping tolerance")
    c.cycles>0 && c.stride>0 || throw(ArgumentError("cycles and stride must be positive"))
    steps=round(Int,c.window/c.dt)
    steps>=1 && abs(steps*c.dt-c.window)<1e-10*c.window || error("window/dt must be an integer")
    isfinite(c.spinup) && c.spinup>=0 && abs(c.spinup/c.dt-round(c.spinup/c.dt))<1e-8 || error("Invalid spinup/dt")
    max(c.mu_forward,c.mu_backward)*c.dt <= 2 || error("Explicit feedback requires mu*dt <= 2 even in isolation")
    c.scheme in (:standard,:cutoff,:diffusive,:voigt,:damping_reversal) || error("Unknown backward scheme")
    c.integrator in (:ifeuler,:ifrk4,:rk4) || error("Unknown integrator")
    c.initialization in (:zero,:projection,:exact,:explicit) || error("Unknown initialization")
    c.schedule in (:continuous,:thirds) || error("Unknown observation schedule")
    c.window_rule in (:halfopen,:legacyclosed) || error("Unknown window rule")
    all(j->j in 1:3,c.components) && length(c.fractions)==3 && all(x->0<=x<=1,c.fractions) || error("Invalid component windows")
    c.scheme==:voigt && c.alpha<=0 && error("Voigt alpha must be positive")
    if c.model==:lorenz
        c.n==3 && c.scheme==:standard && c.integrator==:rk4 || error("Lorenz has three components and uses standard BFN")
        length(c.reference_initial)==length(c.guess_initial)==3 || error("Lorenz states must have length 3")
        all(isfinite,c.reference_initial) && all(isfinite,c.guess_initial) || error("Nonfinite Lorenz state")
        all(x->isfinite(x)&&x>0,(c.sigma,c.rho,c.beta)) || error("Invalid Lorenz coefficients")
    else
        iseven(c.n) && c.n>=16 || error("Use an even grid of at least 16 points")
        0<=c.observation_cutoff<c.n/2 || error("Observation cutoff exceeds grid")
        0<=c.nonlinear_cutoff<c.n/2 || error("Nonlinear cutoff exceeds grid")
        c.profile in (:low,:high,:mixed) || error("Unknown initial profile")
        (c.model==:nse || 0<c.high_mode<c.n/2) || error("Unresolved high mode")
        c.initialization==:explicit && error("Explicit vector initialization is Lorenz-only")
    end
    c.scheme==:damping_reversal && c.model!=:kdv_damped && error("Damping reversal applies to damped KdV")
    c.model==:nse && c.length!=2pi && error("NSE uses a 2pi-periodic square")
    steps
end

function config_dict(c::Config)
    Dict(string(k)=>(getfield(c,k) isa Symbol ? string(getfield(c,k)) : getfield(c,k)) for k in fieldnames(Config))
end
function config_from_dict(d)
    c=Config()
    for (k,v) in pairs(d)
        key=Symbol(k); hasfield(Config,key) || error("Unknown configuration key $k")
        T=fieldtype(Config,key)
        value=T==Symbol ? Symbol(v) : (T<:AbstractVector ? T(collect(v)) : convert(T,v))
        setfield!(c,key,value)
    end
    c
end
