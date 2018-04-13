### The hamiltonian for all calculations
### Time to use multiple dispatch for collecting all stuff

immutable Parameters
    hw::Tuple{Float64,Float64,Float64}
    I::Tuple{Float64,Float64,Float64}
    phi::Tuple{Float64,Float64,Float64}
end


divf(x,y) = y==0 ? 0 : x/y

const sigmax = [0 1;1 0] 
const sigmay = [0 -1im;1im 0]
const sigmaz = [1 0;0 -1] 

δ(m,n) = m==n ? 1 : 0

function H(l_,m_,n_,a_,l,m,n,a,p::Parameters)

    hw,I,phi= p.hw,p.I,p.phi

    sigma = (sigmax[a_,a],sigmay[a_,a],sigmaz[a_,a]) # vector of sigma matricies

    bc = (sqrt(l+1)*δ(l_,l+1), sqrt(m+1)*δ(m_,m+1), sqrt(n+1)*δ(n_,n+1))
    ba = (sqrt(l)*δ(l_,l-1),sqrt(m)*δ(m_,m-1),sqrt(n)*δ(n_,n-1))

    E = δ(l_,l)*δ(m_,m)*δ(n_,n)

    ### Collecting terms
    ### Maybe can be made a bit faster

    s::Complex{Float64} = 0

    s += sum(p.hw.*(l+1/2,m+1/2,n+1/2).*E.*δ(a_,a))
    s += sum(sigma.*I.*(bc.+ba).*(δ(m_,m)*δ(n_,n),δ(l_,l)*δ(n_,n),δ(l_,l)*δ(m_,m)))
    s += sum(sigma.*I.*phi)*E
    return s
end

function FillOperator(H::Function,Nx::Int,Ny::Int,Nz::Int,p::Parameters)
    ### Need to do special care if Is becomes nondiognal!
    #Hm[:,:] = 0
    Hm = Array{Complex}(0)
    I = Array{Int}(0)
    J = Array{Int}(0)
    
    s(l,m,n,a) = 2*(Ny+1)*(Nz+1)*l+2*(Nz+1)*m+2*n+a
    for l in 0:Nx    
        for m in 0:Ny
            for n in 0:Nz
                elements = [
                    (l-1,m,n)
                    (l,m,n)
                    (l+1,m,n)
                    (l,m-1,n)
                    (l,m+1,n)
                    (l,m,n-1)
                    (l,m,n+1)
                ]
                
                for (l_,m_,n_) in elements
                    ### A continue routine
                    if (l_<0) || (m_<0) || (n_<0) || (l_>Nx) || (m_>Ny) || (n_>Nz)
                        continue
                    end
                    for a in 1:2
                        for a_ in 1:2
                            sbra = s(l_,m_,n_,a_)
                            sket = s(l,m,n,a)
                            Hmi = H(l_,m_,n_,a_,l,m,n,a,p)
                            if Hmi!=0
                                push!(I,sbra)
                                push!(J,sket)
                                push!(Hm,Hmi)
                            end
                        end
                    end
                end
            end
        end
    end
    return sparse(I,J,Hm)
end

function FillHamiltonian(Nx::Int,Ny::Int,Nz::Int,p::Parameters)
    return FillOperator(H,Nx,Ny,Nz,p)
end


### A simple operators

OPphix(l_,m_,n_,a_,l,m,n,a,p::Parameters) = (sqrt(l+1)*δ(l_,l+1) + sqrt(l)*δ(l_,l-1)) * δ(m_,m)*δ(n_,n)*δ(a_,a)

OPphiy(l_,m_,n_,a_,l,m,n,a,p::Parameters) = (sqrt(m+1)*δ(m_,m+1) + sqrt(m)*δ(m_,m-1)) * δ(l_,l)*δ(n_,n)*δ(a_,a)

OPphiz(l_,m_,n_,a_,l,m,n,a,p::Parameters) = (sqrt(n+1)*δ(n_,n+1) + sqrt(n)*δ(n_,n-1)) * δ(l_,l)*δ(m_,m)*δ(a_,a)

_p = Parameters((0,0,0),(0,0,0),(0,0,0))

@time const Phix = FillOperator(OPphix,N...,_p)
@time const Phiy = FillOperator(OPphiy,N...,_p)
@time const Phiz = FillOperator(OPphiz,N...,_p)
