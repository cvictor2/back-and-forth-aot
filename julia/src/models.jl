struct SpectralGrid{D,PF,PI}
    kx::Array{Float64,D}
    ky::Array{Float64,D}
    k2::Array{Float64,D}
    radius::Array{Float64,D}
    dx::Array{ComplexF64,D}
    dy::Array{ComplexF64,D}
    invlap::Array{Float64,D}
    keep::BitArray{D}
    obs::BitArray{D}
    linear::Array{ComplexF64,D}
    forcing::Array{ComplexF64,D}
    u0::Array{ComplexF64,D}
    volume::Float64
    forward_plan::PF
    inverse_plan::PI
end

function spectral_grid(c::Config)
    n=c.n; idx=vcat(collect(0:n÷2-1),collect(-n÷2:-1))
    modes=Float64.(idx); k=(2pi/c.length).*modes
    if c.model==:nse
        # MATLAB ndgrid: the x coordinate varies along the first array dimension.
        k[n÷2+1]=pi*n/c.length
        kx=repeat(reshape(k,n,1),1,n); ky=repeat(reshape(k,1,n),n,1)
        radius=sqrt.(reshape(Float64.(idx).^2,n,1) .+ reshape(Float64.(idx).^2,1,n))
        dx=im.*kx;dy=im.*ky;dx[n÷2+1,:].=0;dy[:,n÷2+1].=0
        keep=trues(n,n);lo=ceil(Int,n/3);hi=n-floor(Int,n/3)+1
        keep[lo:hi,:].=false;keep[:,lo:hi].=false;keep[1]=false
        d=matread(c.restart_file);psi=ComplexF64.(d["psi_hat_initial"]);fhat=ComplexF64.(d["f_hat"])
        size(psi)==size(fhat) || error("Restart and forcing shapes differ")
        if size(psi)!=(n,n)
            N=size(psi,1)
            c.allow_restart_resample && size(psi)==(N,N) && iseven(N) && N>n || error("Restart resolution must match n; resampling is explicit and smoke-only")
            take=vcat(1:n÷2,N-n÷2+1:N);scale=(n/N)^2
            psi=psi[take,take].*scale;fhat=fhat[take,take].*scale
        end
        u0=psi;volume=c.length^2
    else
        if c.model in (:kdv_damped,:kdv_viscous)
            # The KdV driver assigns zero to its Nyquist multiplier.
            modes[n÷2+1]=0;k[n÷2+1]=0
        end
        kx=k;ky=zeros(n);radius=abs.(modes);dx=im.*kx;dy=zeros(ComplexF64,n)
        keep=BitArray((radius .<= c.nonlinear_cutoff) .& (radius .> 0));keep[n÷2+1]=false
        x=(0:n-1).*(c.length/n);low=c.low_amplitude.*cos.(2pi.*x./c.length)
        high=c.high_amplitude.*cos.(c.high_mode.*2pi.*x./c.length)
        profile=c.profile==:low ? low : c.profile==:high ? high : low+high
        u0=fft(ComplexF64.(profile));fhat=zeros(ComplexF64,n);volume=c.length
        if c.model in (:kdv_damped,:kdv_viscous)
            fhat=fft(ComplexF64.(c.forcing_amplitude.*exp.(cos.(2pi.*x./c.length))));fhat[1]=0
        end
    end
    k2=kx.^2+ky.^2;invlap=zeros(size(k2));nz=k2.>0;invlap[nz].=-1 ./k2[nz]
    obs=BitArray((radius .< c.observation_cutoff) .& (radius .> 0))
    c.full_observations && (obs .= radius .> 0)
    linear=ComplexF64.(-c.nu.*k2 .- c.gamma)
    c.model==:transport && (linear .+= im*c.transport_speed.*kx)
    c.model in (:kdv_damped,:kdv_viscous) && (linear .+= im*c.dispersion^2 .*kx.^3)
    forcing=c.model==:nse ? invlap.*fhat : fhat
    sample=zeros(ComplexF64,size(k2));pf=plan_fft(sample;flags=FFTW.ESTIMATE);pii=plan_ifft(sample;flags=FFTW.ESTIMATE)
    SpectralGrid(kx,ky,k2,radius,dx,dy,invlap,keep,obs,linear,forcing,u0,volume,pf,pii)
end

function nonlinear(z,g::SpectralGrid,c::Config)
    c.model in (:heat,:transport) && return g.forcing
    masked=z.*g.keep
    if c.model==:nse
        # Velocity is (-psi_y, psi_x), vorticity is Delta*psi.
        u=-real.(g.inverse_plan*(g.dy.*masked));v=real.(g.inverse_plan*(g.dx.*masked))
        uv=g.forward_plan*(u.*v);difference=g.forward_plan*(v.^2-u.^2)
        return @. g.invlap*((g.kx^2-g.ky^2)*uv+g.kx*g.ky*difference)+g.forcing
    end
    u=real.(g.inverse_plan*masked);square=g.forward_plan*(u.^2)
    @. -0.5im*g.kx*square+g.forcing
end

function error_values(e,g::SpectralGrid,c::Config)
    q=abs2.(e)
    c.model==:nse && (q .*= abs2.(g.dx)+abs2.(g.dy))
    hq=g.k2 .*q;scale=sqrt(g.volume)/length(e)
    (scale*sqrt(sum(q)),scale*sqrt(sum(q[g.obs])),scale*sqrt(sum(q[.!g.obs])),
     scale*sqrt(sum(hq)),scale*sqrt(sum(hq[g.obs])),scale*sqrt(sum(hq[.!g.obs])))
end
function legacy_values(e,g::SpectralGrid,c::Config)
    a=max.(abs.(e),1e-14);scale=c.length/c.n
    (scale*norm(a),scale*norm(a[g.obs]),scale*norm(a[.!g.obs]))
end
