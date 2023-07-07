global Ref
global eps_inf
global n_0
global m_eff
global gam

load ref.dat

# Default value :  chi_b = 12.606-1;
# Default value : gam = 1.3e13;
# Default value : n_0 = 1.5e25; # par m^3;
# Default value : m_eff = 0.0738

# Bounds

# Cost function test.
X = [11.606,1.5e25,0.0738,1.3e13];
f(X)
