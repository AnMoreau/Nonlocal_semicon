
clear all;

% Wavelength, nanometers.
lambda=5105;
% Spatial period of the window - discretize means periodize.
d=10*lambda;
% Waist of the incident gaussian beam.
% Smaller than lambda = almost ponctual source.
w=5*lambda;
% Angle of incidence of the beam.
theta=45*pi/180;
% Position of the incident beam with respect to the domain :
% 0.5 means centered on the middle of the domain.
C=0.5;
%----------- Structure
structure
% Image characteristics :
% Number of pixels horizontally
%nx=lambda*10;
nx=floor(d/10);
% Number of points for each layer
ny=floor(hauteur/0.33);
% Number of modes retained for the description of the field
% so that the last mode has an amplitude < 1e-3
nmod=floor(0.83660*d/w);
nmod=0

%----------- From now on, no need to modify anything---------------

hauteur=hauteur/d;
l=lambda/d;
w=w/d;
% !! Omega_p est aussi normalisée !!
w_p=w_p*d;
k0=2*pi/l;

En=zeros(sum(ny),nx+1);
Extn=zeros(sum(ny),nx+1);
Exln=zeros(sum(ny),nx+1);
Eztn=zeros(sum(ny),nx+1);
Ezln=zeros(sum(ny),nx+1);

% Total number of layers
g=length(Type);

% Amplitude of the different modes
X=exp(-w^2*pi^2*[-nmod:nmod].^2).*exp(-2*i*pi*[-nmod:nmod]*C);   % amplitude de l onde incidente
figure(7);
plot(abs(X));
% Scattering matrix corresponding to no interface.
T{1}=[0,1;1,0];

for nm=1:2*nmod+1

% The top medium is assumed to be LOCAL.
  alpha=sqrt(epsilon(Type(1)))*k0*sin(theta)+2*pi*(nm-nmod-1);
  gamma(1)=sqrt(epsilon(Type(1))*k0^2-alpha^2);

  for k=1:g-1

    if (Beta(Type(k))==0)
% Layer scattering matrix for a local, dielectric layer
      Kl(k)=0;
      omega(k)=0;
      t=exp(i*gamma(k)*hauteur(k));
      T{2*k}=[0,t;t,0];
    else
% Layer scattering matrix for a metallic layer.
      Kl(k)=sqrt(alpha^2+((w_p(Type(k))^2)*((1/chi_f(Type(k)))+(1/(1+chi_b(Type(k)))))/Beta(Type(k))^2));
      omega(k)=alpha^2*(1/(1+chi_f(Type(k))+chi_b(Type(k)))-1/(1+chi_b(Type(k))))/Kl(k);
      t=exp(i*gamma(k)*hauteur(k));
      l=exp(-Kl(k)*hauteur(k));
      T{2*k}=[0 0 t 0; 0 0 0 l; t 0 0 0;0 l 0 0];
    end

    gamma(k+1)=sqrt(epsilon(Type(k+1))*k0^2-alpha^2);
% Changing the cut of the square root to achieve perfect stability
    if (imag(gamma(k+1))<0)
      gamma(k+1)=-gamma(k+1);
    end
    if (Beta(Type(k+1))~=0)
      Kl(k+1)=sqrt(alpha^2+((w_p(Type(k+1))^2)*((1/chi_f(Type(k+1)))+(1/(1+chi_b(Type(k+1)))))/Beta(Type(k+1))^2));
      omega(k+1)=alpha^2*(1/(1+chi_f(Type(k+1))+chi_b(Type(k+1)))-1/(1+chi_b(Type(k+1))))/Kl(k+1);
    end

    if ((Beta(Type(k))==0)&&(Beta(Type(k+1))==0))
      b1=gamma(k)/epsilon(Type(k));
      b2=gamma(k+1)/epsilon(Type(k+1));
      T{2*k+1}=[b1-b2,2*b2;2*b1,b2-b1]/(b1+b2);
    else
      if ((Beta(Type(k))==0)&&(Beta(Type(k+1))~=0))
% Interface scattering matrix diel -> metal
	b1=gamma(k)/epsilon(Type(k));
	b2=gamma(k+1)/epsilon(Type(k+1));
	T{2*k+1}=[b1-b2+i*omega(k+1),2*b2,2;2*b1,b2-b1+i*omega(k+1),2; 2*i*omega(k+1)*b1,2*i*omega(k+1)*b2,b1+b2+i*omega(k+1)]/(b1+b2-i*omega(k+1));
      else
% Interface scattering matrix metal -> diel
	b1=gamma(k)/epsilon(Type(k));
	b2=gamma(k+1)/epsilon(Type(k+1));
	T{2*k+1}=[b1-b2+i*omega(k),-2,2*b2; -2*i*omega(k)*b1,b2+b1+i*omega(k),-2*i*omega(k)*b2; 2*b1,-2, b2-b1+i*omega(k)]/(b1+b2-i*omega(k));
      end
    end
  end

% Scattering matrix for the last layer
  t=exp(i*gamma(g)*hauteur(g));
  T{2*g}=[0,t;t,0];

% Once the scattering matrixes have been prepared, now let us combine them.

  H{1}=T{2*g};
  A{1}=T{1};

  for j=1:size(T,2)-2
    A{j+1}= cascade(A{j},T{j+1});
    H{j+1}= cascade(T{size(T,2)-j},H{j});
  end
  r=A{j+1}(1,1);
  disp(r)
% And let us compute the intermediate coefficients from the scattering matrixes

  for j=1:size(T,2)-1
    I{j}= intermediaire(A{j},H{size(T,2)-j});
  end

  h=0;

  I{2*g}=zeros(2,2);

% Computation of the field in the different layers for one mode (plane wave)

  t=1;
  E=zeros(sum(ny),1);
  Exl = zeros(sum(ny),1);
  Ext = zeros(sum(ny),1);
  Ezl = zeros(sum(ny),1);
  Ezt = zeros(sum(ny),1);
  for k=1:g
    for m=1:ny(k)
      h=h+hauteur(k)/ny(k);
      	E(t,1)=I{2*k-1}(1,1)*exp(i*gamma(k)*h)+I{2*k}(2+(Beta(Type(k))~=0),1)*exp(i*gamma(k)*(hauteur(k)-h));
        Ext(t,1)=gamma(k)/epsilon(Type(k))*(I{2*k-1}(1,1)*exp(i*gamma(k)*h)-I{2*k}(2+(Beta(Type(k))~=0),1)*exp(i*gamma(k)*(hauteur(k)-h)));
        Ezt(t,1)=alpha/epsilon(Type(k))*(I{2*k-1}(1,1)*exp(i*gamma(k)*h)+I{2*k}(2+(Beta(Type(k))~=0),1)*exp(i*gamma(k)*(hauteur(k)-h)));
        if (Beta(Type(k))~=0)
          disp(size(I{2*k-1}))
          disp(size(I{2*k}))
          # Je suis pas sûr des conventions de signe.
          #Exl(t,1)=I{2*k-1}(2,1)*exp(-Kl(k)*h)+I{2*k}(4,1)*exp(Kl(k)*(hauteur(k)-h));
          #Ezl(t,1)=Kl(k)/alpha*(I{2*k-1}(2,1)*exp(-Kl(k)*h)-I{2*k}(4,1)*exp(Kl(k)*(hauteur(k)-h)));

          Exl(t,1)=I{2*k-1}(2,1)*exp(Kl(k)*h)+I{2*k}(4,1)*exp(Kl(k)*(hauteur(k)-h));
          #disp(I{2*k-1}(2,1))
          #disp(I{2*k}(4,1))
          disp(h*d)
          disp(exp(-Kl(k)*h))
          disp(exp(Kl(k)*(hauteur(k)-h)))

          Ezl(t,1)=Kl(k)/alpha*(I{2*k-1}(2,1)*exp(Kl(k)*h)-I{2*k}(4,1)*exp(Kl(k)*(hauteur(k)-h)));


        end
	      t=t+1;
      end
      h=0;
  end

% For one mode, the image is invariant horizontally
  E=E*exp(i*alpha*[0:nx]./nx);
  Ext=Ext*exp(i*alpha*[0:nx]./nx);
  Exl=Exl*exp(i*alpha*[0:nx]./nx);
  Ezt=Ezt*exp(i*alpha*[0:nx]./nx);
  Exl=Ezl*exp(i*alpha*[0:nx]./nx);

% All the images have to be combined using the proper amplitude for each plane wave.
  En=En+X(nm)*E;
  Extn = Extn+X(nm)*Ext;
  Exln = Exln+X(nm)*Exl;
  Eztn = Eztn+X(nm)*Ezt;
  Ezln = Ezln+X(nm)*Ezl;

end

V=abs(En);

% You have then to visualize the result. May need some optimization, though.

colormap(jet)

figure(1)
En=abs(En);
imagesc(En)

figure(2)
Extn=abs(Extn);
imagesc(Extn)

figure(3)
A=abs(Exln);
tmp = max(max(A))
B=A+(A==0)*tmp;
tmp = min(min(B));
A = A+(A==0)*tmp;
imagesc(A)

figure(4)
Eztn=abs(Eztn);
imagesc(Eztn)

figure(5)
A=abs(Ezln);
tmp = max(max(A))
B=A+(A==0)*tmp;
tmp = min(min(B));
A = A+(A==0)*tmp;
imagesc(A)
