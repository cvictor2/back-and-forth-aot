const DIAGNOSTIC_COLUMNS=["cycleCoordinate","physicalTime","cycle","direction","step","L2total","L2observed","L2unobserved","H1total","H1observed","H1unobserved"]
const CYCLE_COLUMNS=["cycle";DIAGNOSTIC_COLUMNS[6:11]]
function source_manifest()
    paths=String[]
    for (dir,_,files) in walkdir(joinpath(REPO_ROOT,"julia"))
        for f in files
            (endswith(f,".jl") || f in ("Project.toml","Manifest.toml")) && push!(paths,joinpath(dir,f))
        end
    end
    Dict(relpath(p,REPO_ROOT)=>bytes2hex(sha256(read(p))) for p in sort(paths))
end
function write_csv(path,header,a)
    open(path,"w") do io
        println(io,join(header,','))
        for row in eachrow(a);println(io,join((@sprintf("%.17g",x) for x in row),','));end
    end
end
function write_json(path,d)
    open(path,"w") do io;JSON3.pretty(io,JSON3.write(d));println(io);end
end
"""Create an isolated output folder that survives normal process exit."""
function output_directory(parent,prefix)
    mkpath(parent)
    # mktempdir cleans up at exit unless explicitly disabled. Results are permanent.
    mktempdir(parent;prefix=prefix,cleanup=false)
end
function save_result(r,c)
    mkpath(c.output_root);safe=replace(c.name,r"[^A-Za-z0-9_-]"=>"_")
    folder=output_directory(c.output_root,safe*"_")
    # MAT files contain ordinary arrays/dictionaries, not Julia-serialized objects.
    matwrite(joinpath(folder,"result.mat"),Dict("result"=>r);compress=true)
    write_json(joinpath(folder,"config.json"),config_dict(c))
    write_csv(joinpath(folder,"diagnostics.csv"),DIAGNOSTIC_COLUMNS,r["diagnostics"])
    write_csv(joinpath(folder,"cycle_errors.csv"),CYCLE_COLUMNS,r["cycle_errors"])
    folder
end
function replot_result(path::AbstractString;output=joinpath(dirname(path),"replotted"))
    r=matread(path)["result"];c=config_from_dict(JSON3.read(String(r["config_json"])))
    plot_result(r,c,output)
    output
end
function check_coverage(folder,group,selection=:all)
    rows=readlines(joinpath(REPO_ROOT,"FIGURE_MANIFEST.csv"));expected=String[]
    for line in rows[2:end]
        cells=split(line,',');asset=cells[2]
        cells[3]==string(group)||continue
        if group==:one_d && selection!=:all
            selection==:transport && !(startswith(asset,"IT_")||startswith(asset,"VT_")) && continue
            selection==:burgers && !(startswith(asset,"IB_")||startswith(asset,"VB_")) && continue
            selection==:kdv && !startswith(asset,"bfn_") && continue
            selection==:heat && continue
        end
        push!(expected,asset)
    end
    missing=[name*"."*ext for name in expected for ext in ("png","pdf","svg") if !isfile(joinpath(folder,name*"."*ext))]
    isempty(missing)||error("Missing paper exports: "*join(missing,", "))
    expected
end
