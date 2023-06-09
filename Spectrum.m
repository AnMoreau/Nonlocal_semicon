% This program can compute the reflection and transmission coefficients
% as a function of the wavelength.

%clear all
%clf
%addpath('data/:');

% >>>>>>>>>>>>>>>>>>> Parameters <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Working incidence angle in degrees
% theta=60;
% Polarization - 0 means s or TE; 1 means p or TM.
polarization=1;
% Spectral range in LENGTH UNIT
min=2000;
max=10000;
% Number of points
Npoints=5000;
%_____________________________________________________________________

lambda=linspace(min,max,Npoints);

for k=1:Npoints

	[r(k),t(k)]=coefficient(theta*pi/180,lambda(k));

end

% >>>>>>>>>>>>>>>>>>>   Vizualization   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

figure(1)

plot(lambda,abs(r).^2,'linewidth',2),ylabel('Reflexion'),xlabel('Wavelength'),title('Energy reflection coefficient');
