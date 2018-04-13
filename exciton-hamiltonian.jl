### 5 point calculation

function SparseArray(t::Type,f::Function,indicies::Tuple{Int,Int})
    n,m = indicies

    V = Array{t}(0)
    I = Array{Int}(0)
    J = Array{Int}(0)

    for i in 1:n
        for j in 1:m
            val = f(i,j)
            if val!=0
                push!(I,i)
                push!(J,j)
                push!(V,val)
            end
        end
    end

    return sparse(I,J,V)
end

function d1(i,j,y)
    x(p) = y(i+p) - y(i)

    j==i-2 ? -((x(-1)*x(1)*x(2))/(x(-2)*(x(-2)-x(-1))*(x(-2)-x(1))*(x(-2)-x(2)))) :
        j==i-1 ? (x(-2)*x(1)*x(2))/((x(-2)-x(-1))*x(-1)*(x(-1)-x(1))*(x(-1)-x(2))) :
        j==i ? -(1/x(-1))-1/x(1)-1/x(2)-1/x(-2) :
        j==i+1 ? (x(-2)*x(-1)*x(2))/((x(-2)-x(1))*x(1)*(x(1)-x(-1))*(x(1)-x(2))) :
        j==i+2 ? (x(-2)*x(-1)*x(1))/((x(-2)-x(2))*x(2)*(x(2)-x(-1))*(x(2)-x(1))) :
        0
end

function d2(i,j,y)
    x(p) = y(i+p) - y(i)

    j==i-2 ? (2*(x(1)*x(2)+x(-1)*(x(1)+x(2))))/(x(-2)*(x(-2)-x(-1))*(x(-2)-x(1))*(x(-2)-x(2))) : 
        j==i-1 ? -((2*(x(1)*x(2)+x(-2)*(x(1)+x(2))))/((x(-2)-x(-1))*x(-1)*(x(-1)-x(1))*(x(-1)-x(2)))) :
        j==i ? (2*(x(1)*x(2)+x(-1)*(x(1)+x(2))+x(-2)*(x(-1)+x(1)+x(2))))/(x(-2)*x(-1)*x(1)*x(2)) :
        j==i+1 ? -((2*(x(-1)*x(2)+x(-2)*(x(-1)+x(2))))/((x(-2)-x(1))*x(1)*(x(1)-x(-1))*(x(1)-x(2)))) :
        j==i+2 ? -((2*(x(-1)*x(1)+x(-2)*(x(-1)+x(1))))/((x(-2)-x(2))*x(2)*(x(2)-x(-1))*(x(2)-x(1)))) :
        0
end

### There is a bug if I do put d1,d2 inside a function

function Hamiltonian(α,II,G,p0,L,nn)

    sigmax = [0 1;1 0] 
    sigmay = [0 -1im;1im 0]
    sigmaz = [1 0;0 -1]

    Ix,Iy = II
    Gx,Gy = G
    p0x,p0y = p0

    ### One could make theese nicer
    xx(i) = i<= nn//2 ? -(L+1)^((nn/2-i+1)/(nn/2+1)) + 1 : (L+1)^((i-nn/2)/(nn/2+1)) - 1
    x(i) = i<= nn//2 ? -(L+1)^((nn/2-i+1)/(nn/2+1)) + 1 - xx(nn/2)/2 : (L+1)^((i-nn/2)/(nn/2+1)) - 1 - xx(nn/2+1)/2

    ### Laplacian operator


    D2 = SparseArray(Float64,(i,j)->d2(i,j,x),(nn,nn))
    Lap = Gx*kron(D2,speye(nn)) + Gy*kron(speye(nn),D2)

    ### Potential

    V = Float64[]
    for i in 1:nn
        for j in 1:nn
            x_ = x(i)
            y_ = x(j)
            r = sqrt(x_^2+y_^2)
            push!(V,-α/r)
        end
    end
    V = spdiagm(V)

    ### Weyl coupling

    D1 = SparseArray(Float64,(i,j)->d1(i,j,x),(nn,nn))

    Hm = kron(-1/2*Lap + V,eye(2)) + Ix*kron(1im*kron(D1,speye(nn)) + p0x*speye(nn^2),sigmax) + Iy*kron(1im*kron(speye(nn),D1) + p0y*speye(nn^2),sigmay)
end

