
%%%%%%%%%%%%%%%%%%%  Elementos finitos Bidimensional %%%%%%%%%%%%%%%%%%%%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc
%u = @(x,y)= 100*x  
%f  = @(x, y) 0.0;
%Dirichlet = [1,2,5,6];
Ne=2;
Neq = 2;

P = [0,0,0,0,100,100];

%%%%%%%%%%%%%%%%%%%% Tabele de elementos e nos LG  %%%%%%%%%%%%%%%%%%%%%%%

LG = zeros(4,2);
LG = [1,3; 3,5; 4,6 ; 2,4];


%%%%%%%%%%%%%%%%%%%%%% Vetor equacao EQ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 EQ = [0,0,1,2,0,0];
 
  
%%%%%%%%%%%%%%%%%%%%% Coordenadas dos nós %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

XY = zeros (6,2);
XY = [0,0; 0,0.5; 0.5,0; 0.5,0.5; 1,0; 1,0.5];


  
%%%%%%%%%%%%%%%%%%%%%%%%%%  Motangem da matriz local %%%%%%%%%%%%%%%%%%%%%

Q = zeros(2,2);     %%% Tamanha e elementos da matriz de condutividade Q %%%
               
for i = 1:2
        Q(i,i)= 1;
        
end


Pe = 1;           %%%%%%%%%%%%%%%%%%%%%% Peso %%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = sqrt(3)/3;   %%%%%%%%%%%%%%%%%% Ponto de Gauss %%%%%%%%%%%%%%%%%%%%%%
PG = [-w,-w ;w,-w ;w,w ;-w,w];


K = zeros(2,2);
C = zeros(4,2);
for i= 1:Ne
    C(1,:)= XY(LG(1,i),:);
    C(2,:)= XY(LG(2,i),:);
    C(3,:)= XY(LG(3,i),:);
    C(4,:)= XY(LG(4,i),:);

%    Mesma coisa que o codigo acima faz   
%    for j=1:4
%         C(j,:)= XY(LG(j,i),:);
%    end

    %FOR DOS PONTOS DE GAUSS
    Ke = zeros(4,4);
    for j= 1:4;
        B = calcB(PG(j,1), PG(j,2));
        J = B*C;        
        DetJ = det(J);
        T = inv(J);
        D = T*B ;   
        Ke = Ke + (D'* Q * D * Pe)* DetJ;
    end
    for a = 1 : 4         % percorrendo a matriz local%
        for b = 1 : 4
            n1 = LG(a,i);
            n2 = LG(b,i);
            Globali = EQ(n1);
            Globalj = EQ(n2);
              if Globali > 0 && Globalj > 0
                  K(Globali,Globalj) = Ke(a,b) + K(Globali,Globalj);
              end
        end      
    end
end



%%%%%%%%%%%%%%%%%%%%%%%% Vetor forca local %%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = zeros(Neq,1);

for i = 1: Ne
    C(1,:)= XY(LG(1,i),:);
    C(2,:)= XY(LG(2,i),:);
    C(3,:)= XY(LG(3,i),:);
    C(4,:)= XY(LG(4,i),:);

    Flocal = zeros(4,1);
    Ke = zeros(4,4);                %FOR DOS PONTOS DE GAUSS
    for j= 1:4;
        B = calcB(PG(j,1), PG(j,2));
        J = B*C;
        DetJ = det(J);
        T = inv(J);
        D = T*B ;   
        Ke = Ke + (D'* Q * D * Pe)* DetJ;
    end
    
    Plocal = zeros(4,1);
    for j = 1:4
        no = LG(j,i);
        Plocal(j) = P(no);
    end
    

    Flocal = Flocal - Ke * Plocal;
    
    for l = 1:4
        no = LG(l,i);
        Globall=  EQ(no);
        
        if Globall > 0
            
            F(Globall) =  F(Globall) + Flocal(l);
            
        end
    end
    

end



C = K\F

%%%%%%%%%%%%%%%%%%%%%%%    Grafico  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  






























