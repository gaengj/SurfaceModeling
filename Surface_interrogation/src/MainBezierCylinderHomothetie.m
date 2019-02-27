%close all, clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul et visualisation de surfaces de Bezier
%
% Calcul d'un cylindre avec np patches
% Application d'une homothetie [10,1,0.5] a tous les p.c.
%     cylindre ==> ellipse
%
% Les patches sont individuellement evalues en un 
% nombre fixe de parametres (u,v).
%
% Parametres: 
% num_p : nombre de valeurs de parametres = nombre de points de controle
% num_n : nombre de normales calcules (pour le calcul des isophotes)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_p=20;                    % nombre de valeurs de parametre en direction u et v
num_n=100;		     % plus num_p est grand plus la surface paraitra lisse
			     % et plus le calcul sera long

% Matrice B des points de controle
np=8 % segments on cicle
alpha=linspace(0,2*pi,np+1)

for i=1:np
PA = 1*cos(alpha(i));
PB = 1*cos(alpha(i+1));
PC = 4*cos(alpha(i));
PD = 4*cos(alpha(i+1));
B(1,1,1,i)= PC;
B(1,4,1,i)= PA;
B(1,2,1,i)= (2*B(1,1,1,i)+  B(1,4,1,i))/3;
B(1,3,1,i)= (  B(1,1,1,i)+2*B(1,4,1,i))/3;

% B(2,1,1,i)= (2*PC+PD)/3;
% B(2,4,1,i)= (2*PA+PB)/3;
B(2,1,1,i)= 4*cos((2*alpha(i)+  alpha(i+1))/3);
B(2,4,1,i)=   cos((2*alpha(i)+  alpha(i+1))/3);
B(2,2,1,i)= (2*B(2,1,1,i)+  B(2,4,1,i))/3;
B(2,3,1,i)= (  B(2,1,1,i)+2*B(2,4,1,i))/3;

% B(3,1,1,i)= (PC+2*PD)/3;
% B(3,4,1,i)= (PA+2*PB)/3;
B(3,1,1,i)= 4*cos((  alpha(i)+2*alpha(i+1))/3);
B(3,4,1,i)=   cos((  alpha(i)+2*alpha(i+1))/3);
B(3,2,1,i)= (2*B(3,1,1,i)+  B(3,4,1,i))/3;
B(3,3,1,i)= (  B(3,1,1,i)+2*B(3,4,1,i))/3;

B(4,1,1,i)= PD;
B(4,4,1,i)= PB;
B(4,2,1,i)= (2*B(4,1,1,i)+  B(4,4,1,i))/3;
B(4,3,1,i)= (  B(4,1,1,i)+2*B(4,4,1,i))/3;


PA = 1*sin(alpha(i));
PB = 1*sin(alpha(i+1));
PC = 4*sin(alpha(i));
PD = 4*sin(alpha(i+1));
B(1,1,2,i)= PC;
B(1,4,2,i)= PA;
B(1,2,2,i)= (2*B(1,1,2,i)+  B(1,4,2,i))/3;
B(1,3,2,i)= (  B(1,1,2,i)+2*B(1,4,2,i))/3;

% B(2,1,2,i)= (2*PC+PD)/3;
% B(2,4,2,i)= (2*PA+PB)/3;
B(2,1,2,i)= 4*sin((2*alpha(i)+  alpha(i+1))/3);
B(2,4,2,i)=   sin((2*alpha(i)+  alpha(i+1))/3);
B(2,2,2,i)= (2*B(2,1,2,i)+  B(2,4,2,i))/3;
B(2,3,2,i)= (  B(2,1,2,i)+2*B(2,4,2,i))/3;

% B(3,1,2,i)= (PC+2*PD)/3;
% B(3,4,2,i)= (PA+2*PB)/3;
B(3,1,2,i)= 4*sin(  (alpha(i)+2*alpha(i+1))/3);
B(3,4,2,i)=   sin(  (alpha(i)+2*alpha(i+1))/3);
B(3,2,2,i)= (2*B(3,1,2,i)+  B(3,4,2,i))/3;
B(3,3,2,i)= (  B(3,1,2,i)+2*B(3,4,2,i))/3;

B(4,1,2,i)= PD;
B(4,4,2,i)= PB;
B(4,2,2,i)= (2*B(4,1,2,i)+  B(4,4,2,i))/3;
B(4,3,2,i)= (  B(4,1,2,i)+2*B(4,4,2,i))/3;

B(1,1,3,i)= 0;
B(1,2,3,i)= 3;
B(1,3,3,i)= 6;
B(1,4,3,i)= 9;

B(2,1,3,i)= 0;
B(2,2,3,i)= 3;
B(2,3,3,i)= 6;
B(2,4,3,i)= 9;

B(3,1,3,i)= 0;
B(3,2,3,i)= 3;
B(3,3,3,i)= 6;
B(3,4,3,i)= 9;

B(4,1,3,i)= 0;
B(4,2,3,i)= 3;
B(4,3,3,i)= 6;
B(4,4,3,i)= 9;
end

% Homotetie
B(:,:,1,:) =   2*B(:,:,1,:);
B(:,:,2,:) = 0.5*B(:,:,2,:);

% La matrice B stocke tous les points de controle de tous les patchs
% B(:,:,:,k) sont tous les points de controle du patch k
% La dimension de B(:,:,:,k) est 4 x 4 x 3, i.e., 16 points de controle
% a 3 coordonnees (x,y,z)

% B(:,:,1,k): x-coordonnes des points de controle du patch k comme matrice 4 x 4
% B(:,:,2,k): y-coordonnes des points de controle du patch k comme matrice 4 x 4
% B(:,:,3,k): z-coordonnes des points de controle du patch k comme matrice 4 x 4

% ------------------------------------
% num_p+1 valeurs de parametres uniformes: 0,1,2,...,num_p en u et v
u = linspace(0,1,num_p); 
v = u; 


%  ------------------------------------
% Cubic Bezier patches 
% for k=1:np
%     S(:,:,:,k)=bezierPatchEval(B(:,:,:,k),u,v); %evaluation du patch k en (num_p+1) x (num_p+1) points
% end


% % ------------------------------------
% % Normal vectors of Cubic Bezier patches 
% u=linspace(0,1,num_n); v=u;  %parametrisation uniforme (num_n+1)x (num_n+1) valeurs de parametre
% for k=1:np
%   N(:,:,:,k)=bezierPatchNormal(B(:,:,:,k),u,v); %vecteurs normal du patch k
% end


% ------------------------------------
% Computing Isophotes
  



% ------------------------------------
% Visualisation d'un patch/surface de Bezier
%  plotBezierPatch3D(B(:,:,:,2),S(:,:,:,2)) % plot d'un seul patch k
  plotBezierSurface3D(B,S)		   % plot de tous les np patches
