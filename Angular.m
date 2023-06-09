% This program can compute the reflection and transmission coefficients
% as a function of the incidence angle, as well as the absorption in a
% specified layer.

%clear all

%addpath('data/:');

% >>>>>>>>>>>>>>>>>>> Parameters <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Number of the layer where the absorption has to be computed
Layer=2;
% Workgin wavelength
lambda=8329;
% Polarization - 0 means s or TE; 1 means p or TM.
polarization=1;
% Angular range in degrees here
min=0;
max=80;
% Number of points
Npoints=200;
%_____________________________________________________________________


% Structure geometrical parameters
structure

Ab=zeros(Npoints,length(Type));
theta=linspace(min,max,Npoints);

for k=1:Npoints
     tmp=theta(k)/180*pi;
	   [r(k),t(k)]=coefficient(tmp,lambda);

end

% >>>>>>>>>>>>>>>>>>>   Vizualization   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

figure(1)

plot(theta,abs(r).^2,'linewidth',2),ylabel('Coefficient'),xlabel('Angle'),title('Energy reflection coefficient');


% For test reasons - R+T+Absorption = 1
% R+T+Ab(:,Layer).'-1
