struct StepKernel{A}
    h::Float64
    E::A
    H::A
    V::Array{Float64}
    gain::Float64
end
function kernel(g,c,h; backward=false,truth=false)
    L=copy(g.linear);V=ones(size(L));gain=truth ? 0.0 : backward ? -c.mu_backward : c.mu_forward
    if backward && !truth
        active=c.filtered ? Float64.(.!g.obs) : ones(size(L))
        if c.scheme==:cutoff
            L .+= c.nu.*g.k2 .*(g.radius.>=c.viscosity_cutoff)
        elseif c.scheme==:diffusive
            if c.model in (:kdv_damped,:kdv_viscous)
                L .+= active.*(c.diffusion_correction.*g.k2 .+ (c.cancel_damping ? c.gamma : 0))
            else
                L .+= active.*((c.viscosity_factor-1).*(-c.nu.*g.k2) .+ (c.damping_factor-1)*(-c.gamma))
            end
        elseif c.scheme==:damping_reversal
            L .+= active.*(2c.gamma)
        elseif c.scheme==:voigt
            V .= 1 ./ (1 .+ c.alpha^2 .* g.k2 .* active);L .*= V
        end
    end
    StepKernel(h,exp.(h.*L),exp.(h/2 .*L),V,gain)
end
function physical_step(z,g,c,k)
    h=k.h;E=k.E;H=k.H;V=k.V
    a=V.*nonlinear(z,g,c)
    c.integrator==:ifeuler && return E.*(z+h.*a)
    b=V.*nonlinear(H.*(z+h/2 .*a),g,c)
    d=V.*nonlinear(H.*z+h/2 .*b,g,c)
    e=V.*nonlinear(E.*z+h.*H.*d,g,c)
    E.*z+h/6 .*(E.*a+2 .*H.*(b+d)+e)
end
function spectral_step(z,u,g,c,k)
    out=physical_step(z,g,c,k)
    feedback=k.h*k.gain .*k.V.*g.obs.*(u-z)
    c.feedback_factor && (feedback .*= k.E)
    out+feedback
end

@inline function lorenz_rhs(v,c)
    x,y,z=v
    (c.sigma*(y-x),-c.sigma*x-y-x*z,-c.beta*z+x*y-c.beta*(c.rho+c.sigma))
end
@inline plus_scaled(v,k,h)=ntuple(i->v[i]+h*k[i],3)
@inline function lorenz_physical(v,c,h)
    a=lorenz_rhs(v,c);b=lorenz_rhs(plus_scaled(v,a,h/2),c)
    d=lorenz_rhs(plus_scaled(v,b,h/2),c);e=lorenz_rhs(plus_scaled(v,d,h),c)
    ntuple(i->v[i]+h/6*(a[i]+2b[i]+2d[i]+e[i]),3)
end
function observation_mask(c,t)
    mask=falses(3)
    for j in c.components
        left=(j-1)*c.window/3;right=left+c.fractions[j]*c.window/3
        mask[j]=c.schedule==:continuous || (c.window_rule==:halfopen ? left<=t<right : ((j==1 ? t>=left : t>left) && t<=right))
    end
    mask
end
function lorenz_errors(v,u,c)
    q=ntuple(i->abs2(v[i]-u[i]),3)
    obs=ntuple(i->i in c.components,3)
    (sqrt(sum(q)),sqrt(sum(q[i] for i in 1:3 if obs[i];init=0.0)),sqrt(sum(q[i] for i in 1:3 if !obs[i];init=0.0)),0.,0.,0.)
end
function run_lorenz_experiment(c,steps)
    truth=zeros(3,steps+1);v=Tuple(c.reference_initial)
    for _ in 1:round(Int,c.spinup/c.dt);v=lorenz_physical(v,c,c.dt);end
    truth[:,1].=v
    for j in 1:steps
        v=lorenz_physical(v,c,c.dt);all(isfinite,v)||error("Reference became nonfinite")
        truth[:,j+1].=v
    end
    initial=copy(truth[:,1]);v=Tuple(c.guess_initial)
    c.initialization==:zero && (v=(0.,0.,0.))
    c.initialization==:exact && (v=Tuple(initial))
    c.initialization==:projection && (v=ntuple(i->i in c.components ? initial[i] : 0.,3))
    guess=collect(v);reconstruction=copy(guess);rows=Vector{Float64}[];ends=Vector{Float64}[];states=Vector{Float64}[]
    push!(rows,[0,0,0,0,0,lorenz_errors(v,initial,c)...]);status="completed";failure="";nextidx=1
    for cycle in 1:c.cycles
        failed=false
        for leg in 1:2
            back=leg==2;h=back ? -c.dt : c.dt;gain=back ? -c.mu_backward : c.mu_forward
            indices=back ? (steps+1:-1:2) : (1:steps)
            for (j,idx) in enumerate(indices)
                nextidx=idx+(back ? -1 : 1)
                t=(idx-1+(back ? c.backward_observation_shift : 0))*c.dt
                mask=observation_mask(c,t);physical=lorenz_physical(v,c,h)
                v=ntuple(i->physical[i]+h*gain*mask[i]*(truth[i,idx]-v[i]),3)
                err=lorenz_errors(v,view(truth,:,nextidx),c)
                bad=!all(isfinite,v)||!all(isfinite,err)||maximum(err)>c.max_norm
                if j%c.stride==0 || j==steps || bad
                    push!(rows,[cycle-1+(leg-1)/2+j/(2steps),(nextidx-1)*c.dt,cycle,back ? -1 : 1,j,err...])
                end
                if bad;status="failed";failure="Cycle $cycle, pass $leg, step $j exceeded the error limit";failed=true;break;end
            end
            failed && break
        end
        failed && break
        reconstruction=collect(v);err=lorenz_errors(v,initial,c)
        push!(ends,[cycle,err...]);c.save_cycle_states && push!(states,collect(v))
        if c.stop_tolerance!==nothing && err[1]<=c.stop_tolerance;status="tolerance_reached";break;end
    end
    result_dict(c,rows,ends,states,initial,guess,reconstruction,collect(v),copy(truth[:,nextidx]),copy(truth[:,end]),status,failure)
end

function result_dict(c,rows,ends,states,initial,guess,reconstruction,final,finalref,terminal,status,failure)
    Dict{String,Any}("config_json"=>JSON3.write(config_dict(c)),"diagnostics"=>reduce(vcat,permutedims.(rows)),
      "cycle_errors"=>isempty(ends) ? zeros(0,7) : reduce(vcat,permutedims.(ends)),
      "cycle_states"=>states,"initial_reference"=>initial,"initial_guess"=>guess,
      "reconstruction"=>reconstruction,"final_state"=>final,"final_reference"=>finalref,
      "terminal_reference"=>terminal,"completed_cycles"=>length(ends),"status"=>status,"failure"=>failure,
      "norm_convention"=>"absolute integral L2; velocity L2/H1 seminorm for NSE; Euclidean for Lorenz",
      "cycle_convention"=>"one forward pass and one backward pass per cycle",
      "runtime"=>Dict("julia"=>string(VERSION),"created_utc"=>string(now(UTC)),"machine"=>Sys.MACHINE))
end

"""Run a single explicitly configured experiment, optionally saving plots and data."""
function run_experiment(c::Config)
    steps=validate(c)
    if c.model==:lorenz
        r=run_lorenz_experiment(c,steps)
    else
        g=spectral_grid(c);needed=16length(g.u0)*(steps+1)
        needed<=c.max_reference_bytes || error("Reference needs $needed bytes, above max_reference_bytes=$(c.max_reference_bytes). Reduce the grid/window or increase the explicit limit.")
        ktruth=kernel(g,c,c.dt;truth=true);u=copy(g.u0)
        for _ in 1:round(Int,c.spinup/c.dt);u=physical_step(u,g,c,ktruth);end
        all(isfinite,u)||error("Spinup became nonfinite")
        truth=Array{ComplexF64}(undef,length(u),steps+1);truth[:,1]=vec(u)
        for j in 1:steps
            u=physical_step(u,g,c,ktruth);all(isfinite,u)||error("Reference became nonfinite at step $j")
            truth[:,j+1]=vec(u)
        end
        initial=reshape(copy(truth[:,1]),size(u))
        v=c.initialization==:zero ? zero(initial) : c.initialization==:projection ? initial.*g.obs : copy(initial)
        guess=copy(v);reconstruction=copy(v);rows=Vector{Float64}[];ends=Vector{Float64}[];states=typeof(v)[]
        push!(rows,[0,0,0,0,0,error_values(v-initial,g,c)...]);status="completed";failure="";nextidx=1
        for cycle in 1:c.cycles
            failed=false
            for leg in 1:2
                back=leg==2;h=back ? -c.dt : c.dt;k=kernel(g,c,h;backward=back)
                indices=back ? (steps+1:-1:2) : (1:steps)
                for (j,idx) in enumerate(indices)
                    nextidx=idx+(back ? -1 : 1);old=v
                    v=spectral_step(v,reshape(view(truth,:,idx),size(u)),g,c,k)
                    err=error_values(v-reshape(view(truth,:,nextidx),size(u)),g,c)
                    bad=!all(isfinite,v)||!all(isfinite,err)||maximum(err)>c.max_norm
                    if c.model in (:kdv_damped,:kdv_viscous)
                        bad |= maximum(abs,g.inverse_plan*v)>c.max_amplitude
                        if bad;v=old;nextidx=idx;end
                    end
                    if bad
                        status="failed";failure="Cycle $cycle, pass $leg, step $j exceeded the finite/amplitude/error limit";failed=true
                    end
                    if (!bad || !(c.model in (:kdv_damped,:kdv_viscous))) && (j%c.stride==0 || j==steps || bad)
                        push!(rows,[cycle-1+(leg-1)/2+j/(2steps),(nextidx-1)*c.dt,cycle,back ? -1 : 1,j,err...])
                    end
                    failed && break
                end
                failed && break
            end
            failed && break
            reconstruction=copy(v);err=error_values(v-initial,g,c)
            push!(ends,[cycle,err...]);c.save_cycle_states && push!(states,copy(v))
            if c.stop_tolerance!==nothing && err[1]<=c.stop_tolerance;status="tolerance_reached";break;end
        end
        r=result_dict(c,rows,ends,states,initial,guess,reconstruction,v,reshape(copy(truth[:,nextidx]),size(u)),reshape(copy(truth[:,end]),size(u)),status,failure)
    end
    r["source_manifest_json"]=JSON3.write(source_manifest())
    if c.save_output
        folder=save_result(r,c);r["output_folder"]=folder
        c.make_plots && plot_result(r,c,folder)
    end
    println(c.name,": ",r["status"],", ",r["completed_cycles"]," completed cycles")
    r
end
