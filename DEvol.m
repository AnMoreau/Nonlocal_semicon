% Differential Evolution -- en version pas (trop) verbeuse, pas (trop) surveillable
% Utilisation :
% DEvol(@f,budget,Xmin,Xmax,population)
%  où @f est la fonction de coût (ne pas oublier le '@')
%  budget, ben, le budget alloué à l'optimisation en nombre d'évaluations de la fonction de coût
%  Xmin et Xmax deux vecteurs de même taille (attention) contenant les bornes min et max du domaine d exploration
%   et qui sont utilisés pour le tirage de la génération initiale.
%  population, la taille de la population (30 de façon standard)

function [best,convergence]=DEvol(f_cout,budget,Xmin,Xmax,tribu)

% Paramètres de DE - paramètres potentiels de la fonction
  cr=0.5; % Chances de passer les paramètres du parent à son rejeton.
  f1=0.9;
  f2=0.8;

% Tirage initial dans les bornes fixées
  X_min=repmat(Xmin,tribu,1);
  X_max=repmat(Xmax,tribu,1);
  X=X_min+(X_max-X_min).*rand(tribu,length(Xmin));
% Evaluation + meilleur individu
  for j=1:tribu
      f(j)=f_cout(X(j,:));
  end
  [fmin,ibest]=min(f);
  best=X(ibest,:);

% Initialisations
  evaluation=tribu
  convergence=zeros(1,round(budget/5));
  generation=1;
  convergence(generation)=fmin;

% Boucle DE
  while evaluation<budget
    for j=1:tribu
      crossover=rand(1,length(Xmin))<cr;
      Y=X(j,:)+f1*(X(randi(tribu),:)-X(randi(tribu),:))+f2*(best-X(j,:));
      Y=crossover.*X(j,:)+(1-crossover).*Y;
      if prod((Y>Xmin).*(Y<Xmax))
        ftmp=f_cout(Y);
        evaluation=evaluation+1;
        if (ftmp<f(j))
      	    f(j)=ftmp;
	          X(j,:)=Y;
        end
      end
    end

    generation=generation+1;
    disp("Generation:")
    disp(generation)
#    disp(evaluation)

    [fmin,ibest]=min(f);
    best=X(ibest,:);
    convergence(generation)=fmin;
    disp(fmin)
    disp(best)
  end

  convergence=convergence(1:generation);
%  X=best;
%  save -text file.txt X convergence
endfunction
