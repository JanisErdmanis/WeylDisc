I = (1,1)
G = (1,2)
αs = linspace(1/6,4,100)

N = 400
h = 1e-6 
σoff = -0.2

NEV = 4
L = 20

@everywhere include("exciton-hamiltonian.jl")

######################### Just a Test #####################

@everywhere function energies(α,I,G,p0,L,N,NEV,σoff)
    BLAS.set_num_threads(1)
    
    σ = σoff - 2*α^2 - maximum(I)^2/2 - maximum(I)*maximum(p0)

    println("σ = $σ")
    println("p0 = $p0")
    run(`hostname`)

    @time Hm = Hamiltonian(α,I,G,p0,L,N)
    @time solution = eigs(Hm,nev=NEV,tol=1e-6,maxiter=10000,which=:LM,sigma=σ)

    return real(solution[1])
end

p0 = [((0,0,0),α) for α in αs]
px = [((h,0,0),α) for α in αs]
py = [((0,h,0),α) for α in αs]

@sync begin
    @async global E0 = pmap(x->energies(x[2],I,G,x[1],L,N,NEV,σoff),p0)
    @async global Ex = pmap(x->energies(x[2],I,G,x[1],L,N,NEV,σoff),px)
    @async global Ey = pmap(x->energies(x[2],I,G,x[1],L,N,NEV,σoff),py)
end

DEx = (Ex - E0)/h
DEy = (Ey - E0)/h

using JLD
save("outdir/derivatives-$N-$NEV.jld","αs",αs,"DEx",DEx,"DEy",DEy)
