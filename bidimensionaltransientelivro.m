
%%%%%%%%%%%%%%%%%%%  Elementos finitos Bidimensional %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
%u = @(x,y,t) sin(pi*x).*sin(pi*y)*exp(t);
%f  = @(x, y) 0.0;
%Dirichlet = [1,2,5,6];
ro = 1;
c  = 1;
Nt =  16;
deltaT = 1/Nt;
alpha = 1/2;
Ne = 16;

P = zeros(25,1);
P(10)= 100;
P(15)= 100;
P(20)= 100;
P(22)= 100;
P(23)= 100;
P(24)= 100;
P(25)= 100;


%%%%%%%%%%%%%%%%%%%% Tabele de elementos e nos LG  %%%%%%%%%%%%%%%%%%%%%%%

% LG = zeros(4,16);
LG = [1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19;
    2,3,4,5,7,8,9,10,12,13,14,15,17,18,19,20;
    7,8,9,10,12,13,14,15,17,18,19,20,22,23,24,25;
    6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24];



%%%%%%%%%%%%%%%%%%%%%% Vetor equacao EQ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EQ = [0,0,0,0,0,0,1,2,3,0,0,4,5,6,0,0,7,8,9,0,0,0,0,0,0];

Neq = max(EQ);

%%%%%%%%%%%%%%%%%%%%% Coordenadas dos nós %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% XY = zeros (25,2);
XY = [0,0;
    0.25,0;
    0.5,0;
    0.75,0;
    1,0;
    0,0.25;
    0.25,0.25;
    0.5,0.25;
    0.75,0.25;
    1,0.25;
    0,0.5;
    0.25,0.5;
    0.5,0.5;
    0.75,0.5;
    1,0.5;
    0,0.75;
    0.25,0.75;
    0.5,0.75;
    0.75,0.75
    1,0.75;
    0,1;
    0.25,1;
    0.5,1;
    0.75,1;
    1,1];


%%%%%%%%%%%%%%%%%%%%%%%%%%  Motangem da matriz local %%%%%%%%%%%%%%%%%%%%%

Q = zeros(2,2);     %%% Tamanha e elementos da matriz de condutividade Q %%%

for i = 1:2
    Q(i,i)= 1;
end



Pe = 1;                 %%%%%%%%%%%%%%%%%%%%%% Peso %%%%%%%%%%%%%%%%%%%%%%%
w = sqrt(3)/3;          %%%%%%%%%%%%%%%%%% Ponto de Gauss %%%%%%%%%%%%%%%%%
PG = [-w,-w ;w,-w ;w,w ;-w,w];


K = zeros(Neq,Neq);
C = zeros(4,2);

for i= 1:Ne
    
    C(1,:)= XY(LG(1,i),:);
    C(2,:)= XY(LG(2,i),:);
    C(3,:)= XY(LG(3,i),:);
    C(4,:)= XY(LG(4,i),:);
    
    
    Ke = zeros(4,4);            %%%%%%%%%%% Pontos de Gauss %%%%%%%%%%%%%%
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
% A K bate com a do Vitor visto em 8/1/18


%%%%%%%%%%%%%%%%%%%%%%%% Vetor forca local %%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = zeros(Neq,1);

for i = 1: Ne
    C(1,:)= XY(LG(1,i),:); % Coordenadas dos pontos dos elementos finitos
    C(2,:)= XY(LG(2,i),:);
    C(3,:)= XY(LG(3,i),:);
    C(4,:)= XY(LG(4,i),:);
    
    Flocal = zeros(4,1);
    Ke = zeros(4,4);   % Pontos de Gauss
    
    for j= 1:4;
        
        B = calcB(PG(j,1), PG(j,2));  % Derivadas nos pontos de Gauss
        J = B*C;                      % Produto das coordenadas pelas derivadas nos pontos de Gauss
        DetJ = det(J);
        T = inv(J);
        D = T*B;
        Ke = Ke + (D'* Q * D * Pe) * DetJ;
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

% A F bate com a do Vitor visto em 8/1/18

%%%%%%%%%%%%%%%%% Matriz de Capacidade %%%%%%%%%%%%%%%%%%%%%%%%%

M = zeros(Neq,Neq); %% M Global %%

for i = 1:Ne
    
    C(1,:)= XY(LG(1,i),:);
    C(2,:)= XY(LG(2,i),:);
    C(3,:)= XY(LG(3,i),:);
    C(4,:)= XY(LG(4,i),:);
    
    Me = zeros(4,4); %%%% M local
    
    for j= 1:4
        
        A = calcA(PG(j,1), PG(j,2));
        B = calcB(PG(j,1), PG(j,2));
        J = B*C;
        DetJ = det(J);
        Me = Me + (A * A') * DetJ * ro * c ;
        
    end
    
    for a = 1 : 4         % percorrendo a matriz local%
        for b = 1 : 4
            n1 = LG(a,i);
            n2 = LG(b,i);
            Globali = EQ(n1);
            Globalj = EQ(n2);
            if Globali > 0 && Globalj > 0
                M(Globali,Globalj) = Me(a,b) + M(Globali,Globalj);
            end
        end
    end
end
% A M bate com a do Vitor visto em 8/1/18


%%%%%%%%%%%%%%%%%% Método Trapezoidal Generalizado %%%%%%%%%%%%%%%%%%%%
nos = [1:length(EQ)]';

nosNaoPresc = nos(EQ~=0);
nosPresc    = nos(EQ==0) ;
matSolNum = zeros(length(EQ),Nt+1);

Dn = P(nosNaoPresc);
matSolNum(nosNaoPresc,1) = Dn;
matSolNum(:,1)    = P;

Vn = M \(F - (K * Dn)); % Calculando a velocidade inicial

% Set up the movie.
writerObj = VideoWriter('teste.avi'); % Name it.
writerObj.FrameRate = 1; % How many frames per second.
open(writerObj);

[X,Y] = meshgrid(0:0.25:1, 0:0.25:1);
i=0;
Z = [matSolNum(1,i+1),matSolNum(2,i+1),matSolNum(3,i+1),matSolNum(4,i+1),matSolNum(5,i+1);
    matSolNum(6,i+1),matSolNum(7,i+1),matSolNum(8,i+1),matSolNum(9,i+1),matSolNum(10,i+1);
    matSolNum(11,i+1),matSolNum(12,i+1),matSolNum(13,i+1),matSolNum(14,i+1),matSolNum(15,i+1);
    matSolNum(16,i+1),matSolNum(17,i+1),matSolNum(18,i+1),matSolNum(19,i+1),matSolNum(20,i+1);
    matSolNum(21,i+1),matSolNum(22,i+1),matSolNum(23,i+1),matSolNum(24,i+1),matSolNum(25,i+1)];

surf(X,Y,Z)
frame = getframe(gcf);
writeVideo(writerObj, frame);

for i = 1:Nt
    Dn_Mais1_pred = Dn + (deltaT*(1 - alpha)* Vn);
    Vn_Mais1 = (M  + (alpha * deltaT * K))\( F - (K * Dn_Mais1_pred));
    Dn_Mais1 = Dn_Mais1_pred + (alpha*deltaT*Vn_Mais1);
    
    Vn = Vn_Mais1;
    Dn = Dn_Mais1;
    
    matSolNum(nosNaoPresc,i+1) = Dn;
    matSolNum(nosPresc,i+1)    = P(nosPresc);
    
    Z = [matSolNum(1,i+1),matSolNum(2,i+1),matSolNum(3,i+1),matSolNum(4,i+1),matSolNum(5,i+1);
        matSolNum(6,i+1),matSolNum(7,i+1),matSolNum(8,i+1),matSolNum(9,i+1),matSolNum(10,i+1);
        matSolNum(11,i+1),matSolNum(12,i+1),matSolNum(13,i+1),matSolNum(14,i+1),matSolNum(15,i+1);
        matSolNum(16,i+1),matSolNum(17,i+1),matSolNum(18,i+1),matSolNum(19,i+1),matSolNum(20,i+1);
        matSolNum(21,i+1),matSolNum(22,i+1),matSolNum(23,i+1),matSolNum(24,i+1),matSolNum(25,i+1)];
    
    surf(X,Y,Z)
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
    
end
close(writerObj); % Saves the movie.

