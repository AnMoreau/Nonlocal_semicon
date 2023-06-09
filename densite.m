function n = densite(lambda)

w_p=2*pi*299792458/(lambda*1e-9)*sqrt(12.606);
epsilon_0 = 8.85418782e-12;
m=0.0738*9.1093837015e-31;
n = w_p^2*m*epsilon_0/(1.602176565e-19)^2;

n_0 = 1.5e25; # par m^3;
w_pth = sqrt(n_0 * (1.602176565e-19)^2/(m*epsilon_0));
lamda_pth = 2*pi*299792458/w_pth*1e9
endfunction
