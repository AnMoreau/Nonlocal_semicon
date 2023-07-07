function cost = f(X)

# The experimental data should be stored in a 2-columns
# array named "Ref", a global variable.
# First column :
# Second column : reflectance of the structure

global Ref
global eps_inf
global n_0
global m_eff
global gam

eps_inf = X(1);
n_0 = X(2);
m_eff = X(3);
gam = X(4);

Npoints = size(Ref)(1)
disp(Npoints)

# Angle of incidence, in degrees.
theta = 45

 r = zeros(Npoints,1);
 t = r;

for k=1:Npoints
	[r(k),t(k)] = coefficient(theta*pi/180,Ref(k,1));
end

cost = norm(Ref(:,2)-abs(r).^2);
