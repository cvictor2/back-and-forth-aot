function suite_lorenz(mode,folder)
    smoke=mode==:smoke;base=preset(:lorenz);base.output_root=folder;base.make_plots=false
    if smoke;base.window=.006;base.cycles=2;base.stride=10;end
    runs=Dict{String,Any}[]
    function group_plot(configs,labels,name)
        rs=Dict{String,Any}[]
        for c in configs;r=run_experiment(c);push!(rs,r);push!(runs,r);end
        error_plot(rs,labels,folder,name)
    end
    configs=Config[];labels=String[]
    for phi in (1.,1e-5,1e-10)
        c=deepcopy(base);c.name="invariant_$phi";c.cycles=smoke ? 2 : 5
        c.reference_initial=[0.,0.,0.];c.guess_initial=[0.,0.,phi];c.mu_forward=c.mu_backward=1000
        push!(configs,c);push!(labels,"φ = $phi")
    end
    group_plot(configs,labels,"lorenz_error")
    for (field,values) in ((:sigma,[10,20,30,40,50,100]),(:rho,[10,20,28,30,40,50,100]),(:beta,[3,6,8,9,12,15,30]))
        configs=Config[];labels=String[]
        for value in values
            c=deepcopy(base);setfield!(c,field,Float64(value));field!=:sigma && (c.sigma=20)
            field==:beta && (c.beta=value/3)
            c.name="$(field)_$value";push!(configs,c);push!(labels,field==:beta ? "b = $value/3" : "$field = $value")
        end
        group_plot(configs,labels,"Lorenz_$field")
    end
    for family in (:gain,:fraction)
        values=family==:gain ? [10,100,1000,10000,100000] : [0,1/8,1/4,1/2,3/4,7/8,15/16,31/32,63/64,1]
        configs=Config[];labels=String[]
        for value in values
            c=deepcopy(base);c.cycles=smoke ? 2 : 10;c.components=[1,2,3];c.schedule=:thirds
            c.mu_forward=c.mu_backward=family==:gain ? value : 1e4
            family==:fraction && (c.fractions=fill(Float64(value),3))
            if mode==:historical;c.window_rule=:legacyclosed;c.backward_observation_shift=-1;end
            c.name="window_$(family)_$value";push!(configs,c);push!(labels,"$family = $value")
        end
        group_plot(configs,labels,family==:gain ? "discrete_data" : "discrete_gaps")
    end
    runs
end
function suite_1d(mode,folder,selection)
    selection in (:all,:transport,:burgers,:heat,:kdv)||error("Unknown 1D selection")
    smoke=mode==:smoke;runs=Dict{String,Any}[]
    specs=selection==:heat ? [(:heat,"HT",false)] : [(:transport,"IT",false),(:transport,"VT",true),(:burgers,"IB",false),(:burgers,"VB",true)]
    for (model,tag,viscous) in specs
        selection in (:all,:heat,model)||continue
        for (profile,suffix) in ((:high,"none"),(:low,"full"),(:mixed,"partial"))
            c=preset(model,profile);c.name=tag*"_"*suffix;c.output_root=folder;c.make_plots=false
            if viscous;c.nu=.05;c.scheme=:cutoff;end
            tag=="VT" && profile==:high && (c.high_amplitude=1)
            if smoke;c.n=256;c.nonlinear_cutoff=42;c.dt=1e-5;c.window=.0002;c.spinup=.0001;c.cycles=2;c.stride=2;end
            r=run_experiment(c);push!(runs,r)
            error_plot([r],[c.name],folder,tag*"_error_"*suffix;columns=[6,7,8],cycle_limit=c.cycles)
            initial_plots(r,c,folder;state_name=tag*"_soln_"*suffix,spectrum_name=tag*"_spec_"*suffix)
        end
    end
    if selection in (:all,:kdv)
        specs=((:kdv_damped,:standard,"bfn_normal"),(:kdv_damped,:damping_reversal,"bfn_damped"),
         (:kdv_damped,:diffusive,"bfn_diffusive"),(:kdv_damped,:voigt,"bfn_VOIGT"),
         (:kdv_viscous,:standard,"bfn_model_diffusive"),(:kdv_viscous,:diffusive,"bfn_viscous_model_diffusive"),
         (:kdv_viscous,:voigt,"bfn_VOIGT_model_diffusive"))
        for (model,scheme,name) in specs
            c=preset(model);c.scheme=scheme;c.name=name;c.output_root=folder;c.make_plots=false
            if smoke;c.n=128;c.nonlinear_cutoff=floor(Int,c.n/3);c.window=.0002;c.cycles=2;c.stride=2;end
            r=run_experiment(c);push!(runs,r)
            error_plot([r],[c.name],folder,name;columns=[6,7,8],cycle_limit=c.cycles)
        end
    end
    runs
end
function suite_nse(mode,folder)
    base=preset(:nse);base.output_root=folder;base.make_plots=false
    if mode==:smoke
        base.n=64;base.nonlinear_cutoff=21;base.observation_cutoff=12;base.window=.004;base.cycles=2;base.allow_restart_resample=true
    end
    names=["BFN","dBFN","VBFN","filtered_dBFN","filtered_VBFN","full_BFN","full_dBFN","full_VBFN"]
    schemes=[:standard,:diffusive,:voigt,:diffusive,:voigt,:standard,:diffusive,:voigt]
    runs=Dict{String,Any}[];ends=Float64[]
    for i in eachindex(names)
        c=deepcopy(base);c.name=names[i];c.scheme=schemes[i];c.filtered=i in (4,5);c.full_observations=i>=6
        r=run_experiment(c);push!(runs,r)
        write_csv(joinpath(r["output_folder"],"cycle_errors.csv"),CYCLE_COLUMNS,r["cycle_errors"])
    end
    nse_paper_plots(runs,base,folder)
    open(joinpath(folder,"nse_endpoints.csv"),"w") do io
        println(io,"method,L2,observedL2,unobservedL2,H1,observedH1,unobservedH1")
        for (name,r) in zip(names,runs)
            values=isempty(r["cycle_errors"]) ? fill(NaN,6) : r["cycle_errors"][end,2:7]
            println(io,name,',',join((@sprintf("%.17g",x) for x in values),','))
        end
    end
    runs
end
"""Run a paper-resolution suite, a short smoke suite, or redraw the saved paper curves."""
function run_suite(group::Symbol;mode::Symbol=:paper,selection::Symbol=:all,output=nothing)
    group in (:lorenz,:one_d,:nse)||error("Unknown suite")
    mode in (:paper,:smoke,:saved,:historical)||error("Mode must be paper, smoke, saved or historical")
    mode==:historical && group!=:lorenz && error("Historical mode is Lorenz-only")
    if mode==:saved
        folder=output===nothing ? joinpath(REPO_ROOT,"results","julia","saved_$(group)") : output
        return replot_paper(;group=group,output=folder)
    end
    root=joinpath(REPO_ROOT,"results","julia");mkpath(root)
    folder=output===nothing ? output_directory(root,"$(group)_$(mode)_") : String(output);mkpath(folder)
    runs=group==:lorenz ? suite_lorenz(mode,folder) : group==:one_d ? suite_1d(mode,folder,selection) : suite_nse(mode,folder)
    write_json(joinpath(folder,"runs.json"),[Dict("name"=>JSON3.read(r["config_json"])["name"],"folder"=>get(r,"output_folder",""),"status"=>r["status"],"completed_cycles"=>r["completed_cycles"]) for r in runs])
    check_coverage(folder,group,selection)
    files=sort([splitext(f)[1] for f in readdir(folder) if endswith(f,".svg")])
    write_json(joinpath(folder,"coverage.json"),Dict("group"=>string(group),"mode"=>string(mode),"assets"=>files,"count"=>length(files)))
    println("Saved suite in ",folder)
    runs
end
function cli(group,args)
    length(args)<=2||error("Usage: runner.jl [paper|smoke|saved|historical] [1D selection]")
    mode=isempty(args) ? :paper : Symbol(args[1]);selection=length(args)>1 ? Symbol(args[2]) : :all
    group!=:one_d && selection!=:all && error("A selection is supported only by run_1d.jl")
    run_suite(group;mode=mode,selection=selection)
end
