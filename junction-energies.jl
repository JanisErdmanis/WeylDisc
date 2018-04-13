AXIS = ["y","z"]
NEV = 4

gx,gy,gz = 3,2,1
αx,αy,αz = 1,1,1
Cx,Cy,Cz = (1,1,1).*100

hw = (sqrt(gx/Cx),sqrt(gy/Cy),sqrt(gz/Cz))
I = (αx/sqrt(2)/(gx*Cx)^(1/4),αy/sqrt(2)/(gy*Cy)^(1/4),αz/sqrt(2)/(gz*Cz)^(1/4))

ϕ = linspace(0,1*5.4,100)

@everywhere N = (27,25,51)
@everywhere include("junction-hamiltonian.jl")

@everywhere function energies(hw,I,phi,N,NEV)
    println("phi = $phi")
    run(`hostname`)
    BLAS.set_num_threads(1)
    p = Parameters(hw,I,phi)
    @time Hm = FillHamiltonian(N...,p)
    @time solution = eigs(Hm,nev=NEV,tol=1e-18,maxiter=1000,which=:SR)

    spinx = []
    spiny = []
    spinz = []
    for i in 1:length(solution[1])
        up=solution[2][:,i][1:2:end]
        down=solution[2][:,i][2:2:end]

        spinxi = real(up'*down + down'*up)
        spinyi = real(-1im*up'*down + 1im*down'*up)
        spinzi = up'*up - down'*down

        push!(spinx,spinxi)
        push!(spiny,spinyi)
        push!(spinz,spinzi)
    end

    phix = []
    phiy = []
    phiz = []
    for i in 1:length(solution[1])
        phixi = real(solution[2][:,i]'*Phix*solution[2][:,i])
        phiyi = real(solution[2][:,i]'*Phiy*solution[2][:,i])
        phizi = real(solution[2][:,i]'*Phiz*solution[2][:,i])

        push!(phix,phixi)
        push!(phiy,phiyi)
        push!(phiz,phizi)
    end

    return Dict(:EVALS=>real(solution[1]),:spinx=>spinx,:spiny=>spiny,:spinz=>spinz,:phix=>phix,:phiy=>phiy,:phiz=>phiz) #real(solution[1]),spin
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
#knots = vcat(knots...)

#knots = [knotsx;knotsy;knotsz;knotsxy;knotsyz;knotszx]
results = pmap(x->energies(hw,I,x,N,NEV),knots)

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
save("outdir/energies-$N $hw $I.jld","knots",knots,"results",results,parlist...)
