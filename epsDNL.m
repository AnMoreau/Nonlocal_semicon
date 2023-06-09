function [epsilon,chi_b,chi_f,w_p,beta] = epsDNL(lambda,chi_b,n0,gam,beta)

global eps_inf,n0,gam,beta_nl


# Drude Non-Local model.

w=2*pi*299792458/(lambda*1e-9);
hbar = 1.05457182e-34;
# chi_b = 12.606-1;
chi_b = eps_inf;
epsilon_0 = 8.85418782e-12;
gam = 1.3e13;
# n_0 = 1.5e25; # par m^3;
# It is not difficult to turn m* into a true variable
# that can be retrieved too...
m=0.0738*9.1093837015e-31;
w_p = sqrt(n_0 * (1.602176565e-19)^2/(m*epsilon_0));
# gam = 1.3e13;
chi_f = -w_p^2 /(w*(w+i*gam));
beta = sqrt(3/5)*hbar/m* (3*pi^2*n_0)^(1/3)*1e9;
epsilon = 1+chi_b+chi_f;

end
