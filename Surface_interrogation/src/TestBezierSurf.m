%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fichier de génération des points des surfaces 2, 3 et 4
% On construit des morceaux de surfaces en se basant sur 
% les fonctions trogonometriques cos et sin principalement.
%
% On combine ensuite ces petites routines pour creer des 
% surfaces C0 composées de plusieurs morceaux de surfaces.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for k=1:4
    for z=1:4
        s = [k tan(k+z) z];
        save('surface2', 's', '-ascii', '-append');
    end    
end  

% for k=1:4
%     for z=1:4
%         s = [k cos(k+z) z];
%         save('surface3', 's', '-ascii', '-append');
%     end    
% end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pour la surface 4, j'ai fait une combinaison d'un patch avec un cos
% et d'un patch avec un sin, puis j'ai manuellement collé les 4 points de
% controle extremes des deux patchs afin de former un raccord C0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


