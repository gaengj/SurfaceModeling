%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul du champs des normales pour chaque patch k
% 
% Input:
%  - matrix B de 16 points de controle de dim 3
%     chaque point de controle a 3 coordonnes (x,y,z)
%     taille de B: 4x4x3
%     B(:,:,k) keme coordonnes des16  points de controle, k=1,2,3
%     B(i,j,:) les 3 coordonnes du point de controle b_ij
%  - u Vecteur de |u|=length(u) valeurs de parametre en u
%  - v Vecteur de |v|=length(v) valeurs de parametre en v

% Output:
%  - matrix N de taille num_n x num_n x 3 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N = bezierPatchNormal(B,u,v)

N=[];

% Matrice de changement de base
M = [
    1   0   0  0;
   -3   3   0  0;
    3  -6   3  0;
   -1   3  -3  1
    ];

MT = M'; % transpose de M
XU = [];
XV=[];

for i = 1:length(u)
    for j = 1:length(v)
        U  = [1 u(i) u(i)^2 u(i)^3];
        DU = [0 1 2*u(i) 3*(u(i)^2)];
        VT = [1 v(j) v(j)^2 v(j)^3]'; 
        DVT = [0 1 2*v(j) 3*(v(j)^2)]';
                
        for k = 1:3 % evaluation en (u,v) pour les 3 dimensions
            XU(i,j,k) = DU * M * B(:,:,k) * MT * VT;
            XV(i,j,k) = U * M * B(:,:,k) * MT * DVT;
        end     
    end    
end

for i = 1:length(u)
    for j = 1:length(v)
       N(i,j,:) = normalize(cross(XU(i,j,:),XV(i,j,:)));   
    end    
end

