%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul des isophotes. On a choisi la premiere methode
% pour les tracer : on evalue la fonction I(u,v) = < N(u,v),L >
% pour un ensemble assez dense de valeurs (u, v) et on fait un plot3
% des points X^k(u,v) pour lesquels |I(u, v) - c| < epsilon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IS = isophote(B,N,u,v,c)

S(:,:,:)=bezierPatchEval(B(:,:,:),u,v); 

%Parametres initiaux
L = [0 0 1];
I = [];
IS = zeros(length(u)*length(v),3);
compteur = 0;

% Calcul de I(u,v)
for i = 1:length(u)
    for j = 1:length(v)
        %Produit scalaire
        I(i,j)= N(i,j,1)*L(1) + N(i,j,2)*L(2) + N(i,j,3)*L(3);
        if abs(I(i,j)-c)<0.01
            compteur = compteur + 1;
            for k=1:3
                IS(compteur,k) = S(i,j,k);
            end  
        end
    end
end
%On transpose IS
IS = IS(1:compteur, :)';
%Affichage de IS
IS


