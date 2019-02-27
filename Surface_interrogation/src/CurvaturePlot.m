
function K = CurvaturePlot(B,N,u,v)

K = zeros(length(u),length(v));
H = [];
G=[];

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
XUU=[];
XVV=[];
XUV=[];

for i = 1:length(u)
    for j = 1:length(v)
        U  = [1 u(i) u(i)^2 u(i)^3];
        DU = [0 1 2*u(i) 3*(u(i)^2)];
        DUU = [0 0 2 6*u(i)];
        VT = [1 v(j) v(j)^2 v(j)^3]'; 
        DVT = [0 1 2*v(j) 3*(v(j)^2)]';
        DVVT = [0 0 2 6*v(i)]';
                
        for k = 1:3 % evaluation en (u,v) pour les 3 dimensions
            XU(i,j,k) = DU * M * B(:,:,k) * MT * VT;
            XV(i,j,k) = U * M * B(:,:,k) * MT * DVT;
            XUU(i,j,k) = DUU * M * B(:,:,k) * MT * VT;
            XVV(i,j,k) = U * M * B(:,:,k) * MT * DVVT;
            XUV(i,j,k) = DU * M * B(:,:,k) * MT * DVT;
        end     
    end    
end

% %-----------------------------------------------
% % Construction de H
for i = 1:length(u)
    for j = 1:length(v)
        H(i,j)= dot(XUU(i,j,:),N(i,j,:))*dot(XVV(i,j,:),N(i,j,:)) - dot(XUV(i,j,:),N(i,j,:))*dot(XUV(i,j,:),N(i,j,:));
    end
end

% %-----------------------------------------------
% % Construction de G
for i = 1:length(u)
    for j = 1:length(v)
        G(i,j)= dot(XU(i,j,:),XU(i,j,:))*dot(XV(i,j,:),XV(i,j,:)) - dot(XU(i,j,:),XV(i,j,:))*dot(XV(i,j,:),XU(i,j,:));
    end
end


% %-----------------------------------------------
% % Calcul de la Courbure de Gauss
for i = 1:length(u)
    for j = 1:length(v)
        K(i,j)= H(i,j)/G(i,j);
    end
end


end
