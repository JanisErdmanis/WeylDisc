AXIS = ["x","y"]
NEV = 5
N = 400
ΔN = 10

ϕ = linspace(0,0.0005,100)
I = (1,1) #(sqrt(2),1,0).*4
G = (1,2) #(1,1,0)
L = 20
α = 1/6
σ = -0.2 - 2*α^2
#@everywhere
#N = 600 ### What is a discretization of the box

@everywhere include("exciton-hamiltonian.jl")
@everywhere function energies(α,I,G,p0,L,N,NEV,σoff)

    σ = σoff - maximum(I)^2/2 - maximum(I)*maximum(p0)

    println("σ = $σ")
    println("p0 = $p0")
    run(`hostname`)
    BLAS.set_num_threads(1)

    function EXTRAPOLATOR(x,y)
        x1,x2 = x[1],x[2]
        y1,y2 = y[1],y[2]

        A = (x1^3*y2 - x2^3*y1)/(x1^3-x2^3)
        B = (y1-y2)/(x1^3-x2^3)

        return A
    end

    Hm = Hamiltonian(α,I,G,p0,L,N[1])
    @time E1 = eigs(Hm,nev=NEV,tol=1e-6,maxiter=10000,which=:LM,sigma=σ)[1]
    E1 = sort(real(E1))

    Hm = Hamiltonian(α,I,G,p0,L,N[2])
    @time E2 = eigs(Hm,nev=NEV,tol=1e-6,maxiter=10000,which=:LM,sigma=σ)[1]
    E2 = sort(real(E2))

    Eres = EXTRAPOLATOR(1./N,(E1,E2))

    return Dict(:EVALS=>Eres,:E1=>E1,:E2=>E2,:N1=>N[1],:N2=>N[2]) #real(solution[1]),spin
end

knotsx = Iterators.product((ϕ,0:0,0:0)...)
knotsy = Iterators.product((0:0,ϕ,0:0)...)
knotsz = Iterators.product((0:0,0:0,ϕ)...)
knotsxy = zip(((ϕ, ϕ, zeros(ϕ)) ./ sqrt(2))...)
knotsyz = zip(((zeros(ϕ), ϕ, ϕ) ./ sqrt(2))...)
knotszx = zip(((ϕ, zeros(ϕ), ϕ) ./ sqrt(2))...)

knots = []
for i in AXIS
    vi = eval(parse("knots$i"))
    push!(knots,collect(vi))
end

knots = vcat([i[:] for i in knots]...)

results = pmap(x->energies(α,I[1:2],G[1:2],x[1:2],L,(N-ΔN,N+ΔN),NEV,σ),knots)

parlist = []
Lright=0
for i in AXIS
    push!(parlist,"knots$i")
    push!(parlist,eval(parse("knots$i")))
    push!(parlist,"results$i")

    Li = length(eval(parse("knots$i")))
    Lright += Li
    Lleft = 1+Lright-Li
    push!(parlist,results[Lleft:Lright])
end

using JLD
save("outdir/columbenergies-$N $ΔN $(I[1:2]) $α.jld","knots",knots,"results",results,parlist...)
