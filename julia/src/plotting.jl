# Load the headless plotting backend only when plots are requested.
function ensure_plotting()
    if !isdefined(@__MODULE__, :CM)
        @eval import CairoMakie
        @eval const CM=CairoMakie
    end
end
function export_plot(fig,folder,name)
    mkpath(folder)
    for ext in ("png","pdf","svg");CM.save(joinpath(folder,name*"."*ext),fig);end
end
function curve_axis(fig,slot,xlabel,ylabel;logx=false,logy=true,title="")
    CM.Axis(fig[slot...];xlabel=xlabel,ylabel=ylabel,title=title,xscale=logx ? log10 : identity,yscale=logy ? log10 : identity)
end
function draw_line(ax,x,y,label;logx=false,logy=true,style=:solid)
    xx=Float64.(vec(x));yy=Float64.(vec(y));length(xx)==length(yy)||error("Curve lengths differ")
    for i in eachindex(xx)
        (!isfinite(xx[i]) || !isfinite(yy[i]) || (logx && xx[i]<=0) || (logy && yy[i]<=0)) && (yy[i]=NaN)
    end
    CM.lines!(ax,xx,yy;label=String(label),linewidth=2,linestyle=style)
end
function _error_plot(results,labels,folder,name;columns=[6],ylabel="L² error",cycle_limit=nothing)
    fig=CM.Figure(size=(1100,650),fontsize=18)
    ax=curve_axis(fig,(1,1),"Completed-cycle coordinate",ylabel)
    for (i,r) in enumerate(results), (j,col) in enumerate(columns)
        d=r["diagnostics"];label=length(columns)==1 ? labels[i] : ["Total","Observed","Unobserved"][j]
        draw_line(ax,d[:,1],d[:,col],label;style=(:solid,:dash,:dot)[min(j,3)])
    end
    log_limits!(ax,[r["diagnostics"][:,col] for r in results for col in columns])
    cycle_limit===nothing || CM.xlims!(ax,0,cycle_limit)
    CM.axislegend(ax;position=:rt,labelsize=14)
    export_plot(fig,folder,name)
end
function shell_power(z,g,c;energy=false)
    q=abs2.(z)
    c.model==:nse && (q .*= abs2.(g.dx)+abs2.(g.dy))
    q .*= g.volume/length(z)^2 * (energy ? 0.5 : 1)
    power=zeros(floor(Int,maximum(g.radius))+1)
    for i in eachindex(q);power[floor(Int,g.radius[i])+1]+=q[i];end
    power
end
function _initial_plots(r,c,folder;state_name="initial_state",spectrum_name="initial_spectrum")
    c.model==:lorenz && return
    g=spectral_grid(c);z=r["initial_reference"];v=r["reconstruction"]
    if c.model==:nse
        power=shell_power(z,g,c;energy=true)
        fig=CM.Figure(size=(1100,650),fontsize=18);ax=curve_axis(fig,(1,1),"Integer radial wave-number shell","Kinetic energy per shell";logx=true)
        draw_line(ax,collect(0:length(power)-1),power,"Reference";logx=true)
        log_limits!(ax,[power]);CM.vlines!(ax,[c.n/3];color=:gray,linestyle=:dash)
        export_plot(fig,folder,spectrum_name);return
    end
    x=(0:c.n-1).*(c.length/c.n);u=real.(g.inverse_plan*z);estimate=real.(g.inverse_plan*v)
    fig=CM.Figure(size=(1100,650),fontsize=18);ax=curve_axis(fig,(1,1),"x","State at physical time zero";logy=false)
    draw_line(ax,x,u,"Reference";logy=false)
    draw_line(ax,x,estimate,"Last completed backward endpoint";logy=false,style=:dash)
    CM.axislegend(ax;labelsize=14);export_plot(fig,folder,state_name)
    pu=shell_power(z,g,c);pv=shell_power(v,g,c)
    fig=CM.Figure(size=(1100,650),fontsize=18);ax=curve_axis(fig,(1,1),"Integer Fourier mode","Squared L² contribution";logx=true)
    draw_line(ax,collect(0:length(pu)-1),pu,"Reference";logx=true)
    draw_line(ax,collect(0:length(pv)-1),pv,"Last completed backward endpoint";logx=true,style=:dash)
    peak=max(maximum(pu),maximum(pv))
    peak>0 ? CM.ylims!(ax,max(peak*1e-18,1e-300),peak*4) : log_limits!(ax,[pu,pv])
    CM.axislegend(ax;labelsize=14)
    CM.vlines!(ax,[c.observation_cutoff-1,c.nonlinear_cutoff];color=:gray,linestyle=:dash)
    export_plot(fig,folder,spectrum_name)
end
function _plot_result(r,c,folder)
    error_plot([r],[c.name],folder,"errors_1";columns=[6,7,8],cycle_limit=c.cycles)
    if c.model==:nse
        error_plot([r],[c.name],folder,"errors_2";columns=[9,10,11],ylabel="H¹ error",cycle_limit=c.cycles)
    end
    initial_plots(r,c,folder)
end
text_label(x)=replace(x isa AbstractString ? String(x) : x isa AbstractArray ? join(vec(x)," ") : string(x), "L^2"=>"L²", "H^1"=>"H¹")
as_records(x)=x isa AbstractDict ? [x] : vec(x)
as_records(x::MAT.MatlabStructArray)=vec(Array{Dict{String,Any}}(x))
"""Redraw all 52 manuscript assets from their saved numerical curves, without simulations."""
function _replot_paper(;group=:all,output=joinpath(REPO_ROOT,"results","julia","saved_curves"))
    group in (:all,:lorenz,:one_d,:nse)||error("Unknown group")
    records=as_records(matread(joinpath(REPO_ROOT,"data","reference_curves.mat"))["records"]);names=String[]
    for record in records
        name=splitext(basename(String(record["file"])))[1]
        which=startswith(lowercase(name),"nse") ? :nse : (startswith(name,"Lorenz") || startswith(name,"lorenz") || startswith(name,"discrete")) ? :lorenz : :one_d
        group in (:all,which)||continue
        axes=as_records(record["axes"]);fig=CM.Figure(size=(1100,length(axes)>1 ? 1100 : 650),fontsize=18)
        for (j,a) in enumerate(axes)
            logx=occursin("_spec_",name)||name=="nse_initial_spectrum";logy=!occursin("_soln_",name)
            ax=curve_axis(fig,(j,1),text_label(a["xlabel"]),text_label(a["ylabel"]);logx=logx,logy=logy)
            lines=as_records(a["lines"])
            for (i,line) in enumerate(lines)
                label=text_label(line["label"]);isempty(label)&&(label="Curve $i")
                style=occursin("unobserved",lowercase(label)) ? :dot : occursin("observed",lowercase(label)) ? :dash : :solid
                draw_line(ax,line["x"],line["y"],label;logx=logx,logy=logy,style=style)
            end
            logy && log_limits!(ax,[vec(line["y"]) for line in lines])
            if !logx && logy
                xmax=maximum(maximum(line["x"]) for line in lines)
                xmax>0 && CM.xlims!(ax,0,xmax)
            end
            length(lines)>1 && CM.axislegend(ax;position=:rt,labelsize=12)
            if logx;CM.vlines!(ax,name=="nse_initial_spectrum" ? [1024/3] : [15,682];color=:gray,linestyle=:dash);end
        end
        export_plot(fig,output,name);push!(names,name)
    end
    write_json(joinpath(output,"coverage.json"),Dict("group"=>string(group),"assets"=>sort(names),"count"=>length(names),"source"=>"data/reference_curves.mat","mode"=>"redraw saved curves, no simulations"))
    names
end
function _nse_paper_plots(results,c,folder)
    initial_plots(results[1],c,folder;spectrum_name="nse_initial_spectrum")
    for (name,indices) in (("nse_H1_histories",1:3),("nse_filtered_H1_histories",1:5))
        fig=CM.Figure(size=(1100,900),fontsize=18)
        for (j,col) in enumerate((10,11))
            ax=curve_axis(fig,(j,1),"Completed-cycle coordinate",j==1 ? "Observed H¹ error" : "Unobserved H¹ error")
            for i in indices
                r=results[i];d=r["diagnostics"];cfg=config_from_dict(JSON3.read(r["config_json"]))
                draw_line(ax,d[:,1],d[:,col],cfg.name)
            end
            log_limits!(ax,[results[i]["diagnostics"][:,col] for i in indices])
            CM.xlims!(ax,0,c.cycles);CM.axislegend(ax;labelsize=13)
        end
        export_plot(fig,folder,name)
    end
end

function log_limits!(ax,series)
    positive=Float64[]
    for values in series, y in values
        isfinite(y) && y>0 && push!(positive,y)
    end
    isempty(positive) && (CM.ylims!(ax,1e-16,1.0);return)
    lo,hi=extrema(positive);a=log10(lo);b=log10(hi);pad=max(0.2,0.04*(b-a))
    CM.ylims!(ax,10.0^max(-307,a-pad),10.0^min(307,b+pad))
end
function error_plot(args...; kwargs...)
    ensure_plotting()
    Base.invokelatest(_error_plot,args...;kwargs...)
end

function initial_plots(args...; kwargs...)
    ensure_plotting()
    Base.invokelatest(_initial_plots,args...;kwargs...)
end

function plot_result(args...; kwargs...)
    ensure_plotting()
    Base.invokelatest(_plot_result,args...;kwargs...)
end

function replot_paper(args...; kwargs...)
    ensure_plotting()
    Base.invokelatest(_replot_paper,args...;kwargs...)
end

function nse_paper_plots(args...; kwargs...)
    ensure_plotting()
    Base.invokelatest(_nse_paper_plots,args...;kwargs...)
end
