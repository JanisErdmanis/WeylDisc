using PyPlot

# TeX labels
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':6})
#rc('text', usetex=True)
#rc('xtick',labelsize=6)

g1=3.
g2=2.
g3=1.
eps=1.e-3

# calculations for plotting the POTATO stats here

function rmin(u)
    "minimal radius from zeros of the determinant of the Hessian of classical energy"
    u1=u[1]
    u2=u[2]
    u3=u[3]
    s2=g1 * u1^2 + g2 * u2^2 + g3 * u3^2;
    s1=-(g1*g2+g1*g3+g2*g3) + g1*g2*u3^2 + g1*g3*u2^2 + g2*g3*u1^2;
    s0=g1*g2*g3;
    discr=s1^2-4*s0*s2
    if discr > 0
        q0 = (-s1+(s1^2-4*s0*s2)^0.5)/2/s2
        q1 = (-s1-(s1^2-4*s0*s2)^0.5)/2/s2
        qq=min(q0,q1)
        ret=1/qq
        if  ret > 1/g2-eps && ret < 1/g3+eps
            return ret;
        end
    end 
end 

function phi1(X)
    theta=X[1]
    phi=X[2]
    u1=sin(theta)*cos(phi)
    u2=sin(theta)*sin(phi)
    u3=cos(theta)
    r=rmin([u1,u2,u3])
    return (r-1/g1)*u1;
end 

function phi2(X)
    theta=X[1]
    phi=X[2]
    u1=sin(theta)*cos(phi)
    u2=sin(theta)*sin(phi)
    u3=cos(theta)
    r=rmin([u1,u2,u3])
    return (r-1/g2)*u2;
end 

function phi3(X)
    theta=X[1]
    phi=X[2]
    u1=sin(theta)*cos(phi)
    u2=sin(theta)*sin(phi)
    u3=cos(theta)
    r=rmin([u1,u2,u3])
    return (r-1/g3)*u3;
end 

# save data for potato in an array
NN=100
ths=linspace(0,pi/2,NN)
phs=linspace(0,2*pi,NN)

pp=Iterators.product(ths,phs)
P1=map(phi1,pp)

pp=Iterators.product(ths,phs)
P2=map(phi2,pp)

pp=Iterators.product(ths,phs)
P3=map(phi3,pp)

P1=reshape(P1,(NN,NN))
P2=reshape(P2,(NN,NN))
P3=reshape(P3,(NN,NN))

# calculate disc region as a surface
ph1d(X)=(1/g3-1/g1)*X[1]*cos(X[2]); 
ph2d(X)=(1/g3-1/g2)*X[1]*sin(X[2]);
ph3d(X)=0;

rs=linspace(0,1,NN)
ps=linspace(0,2*pi,NN)

pp=Iterators.product(rs,ps)
P1a=map(ph1d,pp)

pp=Iterators.product(rs,ps)
P2a=map(ph2d,pp)

pp=Iterators.product(rs,ps)
P3a=map(ph3d,pp)

P1a=reshape(P1a,(NN,NN))
P2a=reshape(P2a,(NN,NN))
P3a=reshape(P3a,(NN,NN))

# Plotting starts here

fig2 = figure()
fig2[:set_size_inches]([3.487,2*3.487/1.618])
#plt.tick_params(labelsize=6)

#ax = fig2[:add_subplot](111),# projection="3d")
#ax[:set_xlabel]("phi_1^r")
#ax[:set_ylabel]("phi_2^r")
#ax[:set_zlabel]("phi_3^r")

# plot disc region
#plot_surface(P1a,P2a,P3a,alpha=1,color="r")
#ax[:plot_surface](P1a, P2a, P3a,alpha=1,color="r")
# plot potato
#ax.plot_wireframe(P1, P2, P3)
#ax.plot_wireframe(P1l, P2l, P3l)

#ax.plot_surface(P1l, P2l, P3l,rstride=1,cstride=1,color="w")
plot_surface(P1, P2, P3,rstride=4,cstride=2,color="w",alpha=0.8)

ax = gca()
ax[:set_xlim](1/g1-1/g3,1/g3-1/g1)
ax[:set_ylim](1/g2-1/g3,1/g3-1/g2)
ax[:set_zlim](1/g2-1/g3,1/g3-1/g2)


#ax[:view_init](100,30)

savefig("fig2bot.svg")
