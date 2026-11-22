# Optional integration checks. Run after the fast numerical tests.
using BFNPaper, Test, JSON3
mktempdir() do folder
    names=replot_paper(output=joinpath(folder,"saved"))
    @test length(names)==52
    for name in names, ext in ("png","pdf","svg")
        @test filesize(joinpath(folder,"saved",name*"."*ext))>0
    end
    for (group,count,assets) in ((:lorenz,38,6),(:one_d,19,43),(:nse,8,3))
        target=joinpath(folder,string(group))
        runs=run_suite(group;mode=:smoke,output=target)
        @test length(runs)==count
        @test all(r->r["status"]=="completed",runs)
        @test length(BFNPaper.check_coverage(target,group))==assets
    end
end
