%%%%%%%%%%%%%%              Determinar a solucao       %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%     - alfa u_xx(x)+ beta u(x) = f(x)    %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%        Metodo dos Elementos Finitos     %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%             Fluxo - Neumann             %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc

%%%%%%%%%%%%% Entradas %%%%%%%%%%%%%%%%%
Ne   = 20;            % Numero de intervalos  
h    = 1/Ne;         % Tamanho dos intervalos 
alfa = 1;            % Constantes positiva    
beta = 1;            % Constantes positiva    
f    =@(x) x;        % Funcao                 
u    =@(x) x - (exp(x))+ ((exp(1)-1)*(exp(-x)+exp(x)))/(exp(1)-exp(-1));
X = [0:h:1]';
he =h;
p=0;                 
q=0;                

%%%%%%%%%%% Construcao da matrizes rigides local %%%%%%%%%%%%%%%%%%%%%%
Ke = zeros(2,2);                       % Tamanho da matriz Ke %

for i = 1:2 
    Ke(1,1)= alfa/he + beta*he/3;      % Entrada da matriz Ke %
    Ke(1,2)= -alfa/he + beta*he/6;    
    Ke(2,1)= -alfa/he + beta*he/6;
    Ke(2,2)= alfa/he + beta*he/3;
end
    
%%%%%%%%% Construcao da matriz rigida global %%%%%%%%%%%%%%%%%%

K = zeros(Ne+1,Ne+1);                   % Tamanho da matriz K %

for e = 1:Ne
    K(e,e)= K(e,e) + Ke(1,1);           % Entradas da matriz K %
    K(e,e+1)= K(e,e+1) + Ke(1,2);
    K(e+1,e)= K(e+1,e) + Ke(2,1);
    K(e+1,e+1)= K(e+1,e+1) + Ke(2,2);                         
end


%%%%%%%%%%% Construcao do vetor de forca local e global %%%%%%%%%%%%%

Fe = zeros(2,1);                              % Tamanho da Matriz F %
F= zeros(Ne+1,1);                             % Tamanho da Matriz F % 


for e=1:Ne
    Fe(1)= he/6*(2*f(X(e))+f(X(e+1)));         % Entrada da matriz Fe %
    Fe(2)= he/6*(f(X(e))+ 2*f(X(e+1)));    
    F(e)= F(e) + Fe(1);                       % Tamanho da Matriz F % 
    F(e+1)= F(e+1) + Fe(2);
end

F(1)=(he/6)*(2*f(X(1))+f(X(2)))-alfa*p;   
F(Ne+1)=he/6*(f(X(Ne))+ 2*f(X(Ne+1)))+alfa*q; 


%%%%%%%%%%%%%%%%% Encontrando a solucao numerica %%%%%%%%%%%%%

C=K\F;

C

%%%%%%%%%%%%%%%%%%%%%%%    Grafico  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot([0:0.001:1],u([0:0.001:1]),'r')
hold on
plot(X,[C],'*b')
grid on

%e = sqrt(h*(u(X(2:end-1))-C)'*(u(X(2:end-1))-C))

