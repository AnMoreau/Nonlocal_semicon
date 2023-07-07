global Ref
global eps_inf
global n_0
global m_eff
global gam

load ref.dat

# Default value X(1) :  chi_b = 12.606-1;
# Default value X(2) : n_0 = 1.5e25; # par m^3;
# Default value X(3) : m_eff = 0.0738
# Default value X(4) : gam = 1.3e13;
#X = [11.606,1.5e25,0.0738,1.3e13];
# Cost function test.
#f(X)


# Bounds
Xmin = [10,1e24,0.01,1e10]
Xmax = [15,1e26,1.,1e14]

[best,convergence]=DEvol(@f,10000,Xmin,Xmax,30);
