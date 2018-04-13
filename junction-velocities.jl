q = linspace(0.01,5,100)
x = 4*q.^2

g = (3,2,1)
α = (1,1,1)
C = (1,1,1)

@everywhere N = (25,25,50)
h = 0.000001

@everywhere include("junction-hamiltonian.jl")

######################### Just a Test #####################

@everywhere function energies(hw,I,phi,N)
    println("hw = $hw")
    println("I = $I")


    BLAS.set_num_threads(1)
    p = Parameters(hw,I,phi)
    @time Hm = FillHamiltonian(N...,p)
    @time solution = eigs(Hm,nev=4,tol=1e-18,maxiter=1000,which=:SR)
    return real(solution[1])
end


### Making a list of parameters
s = map(x->((sqrt(3/x),sqrt(2/x),sqrt(1/x)),(1/sqrt(2)/(3*x)^(1/4),1/sqrt(2)/(2*x)^(1/4),1/sqrt(2)/(1*x)^(1/4))),x)

@sync begin
     @async global E0 = pmap(s->energies(s...,(0,0,0),N),s)
     @async global Ex = pmap(s->energies(s...,(h,0,0),N),s)
     @async global Ey = pmap(s->energies(s...,(0,h,0),N),s)
     @async global Ez = pmap(s->energies(s...,(0,0,h),N),s)
end

DEx = (Ex - E0)/h
DEy = (Ey - E0)/h
DEz = (Ez - E0)/h

using JLD
save("outdir/derivatives-$N-$q-$g.jld","q",q,"DEx",DEx,"DEy",DEy,"DEz",DEz)
