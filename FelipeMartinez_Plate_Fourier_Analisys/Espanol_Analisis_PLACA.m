% ´◔‿ゝ◔`)━☞ ANÁLISIS DE UNA PLACA CARGADA EN EL CENTRO

% La placa esta apoyada en todos los bordes con apoyos simples

%{
************************* DESCRIPCIÓN ***************************
Una placa de dimensiones "a" por "b", está apoyada en todos sus bordes y
tiene una carga de superfice con forma cuadrada de dimensión "c" por "c",
la carga esta en newtons y las dimensiones en milimetros.

TU DEFINES
1) la carga
2) las dimensiones
3) las orientaciones
4) el material

OBTIENES
1) Gráficas 3D de las curvaturas y deformaciones z sobre el espesor w(x,y)
2) Valores de esfuerzos, deformaciones y curvaturas en un punto de análisis
   que tu elijas y en los bordes de forma predeterminada
3) Factores de resistencia por Tsai-Hill y Tsai-Wu

    Este es un ejercicio donde puedes modificar las entradas como la
    cantidad de capas del material compuesto, correr el programa y ver si
    falla o no el compuesto

                        DIVIERTETE!!!

*******************************************************************
 %}
clear;close all; clear all; clc;

orientaciones = [45 90 -45 30 0 -30 -45 90 45 0];

% Dimensiones de la placa
largo = 160; % mm
ancho = 160;  % mm

capas = length(orientaciones);

%__________________________________________________________________________
%          PROPIEDADES MECÁNICAS DEL MATERIAL (l,t) 

% Glass VE al 60% de fibra 

         El = 31161   ;% MPA
         Et = 7452   ; % MPA
        Glt = 2738   ; % MPA
        vlt = 0.3397   ;
        vtl = vlt*(Et/El);

% RESISTENCIA Máximos Esfuerzos
        slmaxT = 500   ;%   MPA
        slmaxC = 500   ;%   MPA
        stmaxT = 50   ; %   MPA
        stmaxC = 50   ; %   MPA
        taumax = 50   ; %   MPA
        
        t=0.2; % espesor de cada capa [mm]
     


% 1.- Hallar [S]lt  y   [Q]lt
s11=1/El;   s12=-vtl/Et;    s21=-vlt/El;    s22=1/Et;   s66=1/Glt;

Slt=[s11    s12     0;
     s21    s22     0;
     0      0       s66];
 
Qlt = inv(Slt);

q11 = Qlt(1,1);     q12 = Qlt(1,2);     q13 = Qlt(1,3);
q21 = Qlt(2,1);     q22 = Qlt(2,2);     q23 = Qlt(2,3);
q31 = Qlt(3,1);     q32 = Qlt(3,2);     q33 = Qlt(3,3);

% -------------------------------------------------------------------
% 2.- Obtener Sxy y Qxy de cada capa
% 3.- Ciclos for de las capas 
% 4.- Obtener ABD y abd

N = length(orientaciones); % Cantidad de capas
h = zeros (N+1);           % Matriz cuadrada de ceros de (capas+1)x(capas+1)

h(N+1) = N*t/2;            % Cambia el primer valor de la ultima fila
qbarK = zeros(3*N,3);        % Matriz de ceros de tamaño (3*capas)x3 

for k = 1 : N
    c = cosd(orientaciones(k));
    s = sind(orientaciones(k));
    
    % Transformación de la matriz de rigidez reducida
    qactual = zeros(3,3); % Matriz vacia de 3x3 que se va llenando / se limpia
    
    qactual(1,1) = q11*c^4+q22*s^4+s*c*(s*c*(q21+q12+4*q33)+c^2*(2*q31+2*q13)+s^2*(2*q32+2*q23));
    qactual(1,2) = q12*c^4+q21*s^4+s*c*(s*c*(q11+q22-4*q33)+c^2*(2*q32-2*q13)+s^2*(2*q31-2*q23));
    qactual(1,3) = s*c*(s*c*(-2*q31+2*q32)+c^2*(-q11+q12)+s^2*(-q21+q22))+(c^2-s^2)*(q13*c^2+q23*s^2+2*q33*s*c);
    qactual(2,1) = q21*c^4+q12*s^4+s*c*(s*c*(q11+q22-4*q33)+c^2*(-2*q31+2*q23)+s^2*(-2*q32+2*q13));
    qactual(2,2) = q22*c^4+q11*s^4+s*c*(s*c*(q21+q12+4*q33)+c^2*(-2*q32-2*q23)+s^2*(-2*q31-2*q13));
    qactual(2,3) = s*c*(s*c*(2*q31-2*q32)+c^2*(-q21+q22)+s^2*(-q11+q12))+(c^2-s^2)*(q13*s^2+q23*c^2-2*q33*s*c);
    qactual(3,1) = s*c*(s*c*(-2*q13+2*q23)+c^2*(-q11+q21)+s^2*(-q12+q22)+2*q33*(c^2-s^2))+(c^2-s^2)*(q31*c^2+q32*s^2);
    qactual(3,2) = s*c*(s*c*(2*q13-2*q23)+c^2*(-q12+q22)+s^2*(-q11+q21)-2*q33*(c^2-s^2))+(c^2-s^2)*(q31*s^2+q32*c^2);
    qactual(3,3) = s*c*(s*c*(q11-q21-q12+q22)+(c^2-s^2)*(-q31+q32-q13+q23))+q33*(c^2-s^2)^2;
   
    qbarK([3*k-2:3*k],[1:3])=qactual([1:3],[1:3]); % GUARDA TODAS LAS MATRICES qactual's
    h(k)=(k-N/2-1)*t;% Todos los valores de las cotas

end
A = zeros(3);   B = zeros(3);   D = zeros(3);

% qbarK

for i=1:3 % Ciclo que hace la suma de los productos de (Qij)(zk-zk1)²³
    for j=1:3
            qactual([1:3],[1:3])=qbarK([1:3],[1:3]);
            A(i,j) = qactual(i,j) * (h(2) - h(1));
            B(i,j) = 1/2*(qactual(i,j) * (h(2)^2 - h(1)^2));
            D(i,j) = 1/3*(qactual(i,j) * (h(2)^3 - h(1)^3));
            
         for k = 2 : N
            qactual([1:3],[1:3]) = qbarK( [3*k-2:3*k] , [1:3] );
            A(i,j) = qactual(i,j) * (h(k+1) - h(k)) + A(i,j);
            B(i,j) = 1/2*(qactual(i,j) * (h(k+1)^2 - h(k)^2)) + B(i,j);
            D(i,j) = 1/3*(qactual(i,j) * (h(k+1)^3 - h(k)^3)) + D(i,j);    
         end
    end
end

    if rem(length(orientaciones),2)==0
    A= A.*[1 1 0 ; 1 1 0 ; 0 0 1];
    B= B.*[0 0 0 ; 0 0 0 ; 0 0 0];
    D= D.*[1 1 1 ; 1 1 1 ; 1 1 1];
    end
    A
    B
    D

ABD = [A,B;B,D];                    abd = inv(ABD);

% Obtener los elementos Dij
D11 = D(1,1);   D12 = D(1,2);   D13 = D(1,3);
D21 = D(2,1);   D22 = D(2,2);   D23 = D(2,3);
D31 = D(3,1);   D32 = D(3,2);   D33 = D(3,3);

% -------------------------------------------------------------------
% 5.- Sacar w(x,y), kx, ky, kxy
% Identico AL CASO 1, toda la carga distribuida en toda la superficie
   

%CARGA

p0=0.002;% MPa

 % ~~~~~~~~~~~ ENCONTRANDO w(x,y), kx, ky, kxy ~~~~~~~~~~~~
      
 % COEFICIENTES COMUNES
    a=largo;    b=ancho;    R=a/b;  
    
 % NUMERO DE ITERACIONES
    mn=input('Ingrese el número de iteraciones: ');
    m=2*mn;   % m = filas     (i)
    n=m;    % n = columnas  (j)
    
                % Matrices iniciales
    mwxy_c=[];      mwxy_00=[];     mwxy_a0=[];     mwxy_0b=[];     mwxy_ab=[]; 
    mkx=[];     mky=[];    
    mkxy_00=[];     mkxy_a0=[];     mkxy_0b=[];     mkxy_ab=[];     
      
    %   Para kx, ky y w(x,y) en el CENTRO
    
for i = 1: 2: m %Coordenadas (x,y) = (a/2 ,b/2) con paso 2 para i impares
    i;    vwxy=0;    vkx = 0;  vky=0;
    x=a/2;    y=b/2;
    for j = 1: 2 : n
        j;
        vwxy = sin(i*pi()*x/a)*sin(j*pi()*y/b) / (i*j*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
        vkx = i*sin(i*pi()*x/a)*sin(j*pi()*y/b) / (j*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
        vky = j*sin(i*pi()*x/a)*sin(j*pi()*y/b) / (i*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
        mwxy_c(i,j)=vwxy;     mkx(i,j) = vkx;       mky(i,j) = vky; % Insertar valores en la matriz
        mwxy_c;       mkx;        mky;
    end
    vwxy;   vkx;    vky;
end
mwxy_c;   mwxy_c=sum(mwxy_c,'all');      wxy_c=mwxy_c*(16*p0*a^4)/(pi()^6)
mkx;    mkx=sum(mkx,'all');              kx=mkx*(16*p0*a^2)/(pi()^4)
mky;    mky=sum(mky,'all');              ky=mky*(16*p0*a^2*R^2)/(pi()^4)

%   Para w(x,y) y kxy en LOS BORDES

for i = 1 : m % Coordenadas (x,y) = (0 , 0)
    vkxy=0; vwxy=0;
    x=0;    y=0;
    for j = 1 : n
        vwxy = sin(i*pi()*x/a)*sin(j*pi()*y/b) / (i*j*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
        vkxy = cos(i*pi()*x/a)*cos(j*pi()*y/b) / (D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*(R^2*i^2*j^2));
        mkxy_00(i,j) = vkxy;        mwxy_00(i,j)=vwxy; % Insertar valores en la matriz
        mkxy_00;
    end
    vkxy;
end
mkxy_00;   mkxy_00=sum(mkxy_00,'all');      kxy_00=mkxy_00*(-32*R*a^2*p0)/(pi()^4)
mwxy_00;   mwxy_00=sum(mwxy_00,'all');      wxy_00=mwxy_00*(16*p0*a^4)/(pi()^6)


for i = 1 : m % Coordenadas (x,y) = (a , 0)
    vkxy=0;     vwxy=0;
    x=a;    y=0;
    for j = 1 : n
        vwxy = sin(i*pi()*x/a)*sin(j*pi()*y/b) / (i*j*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
        vkxy = cos(i*pi()*x/a)*cos(j*pi()*y/b) / (D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*(R^2*i^2*j^2));
        mkxy_a0(i,j) = vkxy;    mwxy_a0(i,j)=vwxy;% Insertar valores en la matriz
        mkxy_a0;
    end
    vkxy;
end
mkxy_a0;   mkxy_a0=sum(mkxy_a0,'all');      kxy_a0=mkxy_a0*(-32*R*a^2*p0)/(pi()^4)
mwxy_a0;   mwxy_a0=sum(mwxy_a0,'all');       wxy_a0=mwxy_a0*(16*p0*a^4)/(pi()^6)


for i = 1 : m % Coordenadas (x,y) = (0 , b)
    vkxy=0;     vwxy=0;
    x=0;    y=b;
    for j = 1 : n
        vwxy = sin(i*pi()*x/a)*sin(j*pi()*y/b) / (i*j*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
        vkxy = cos(i*pi()*x/a)*cos(j*pi()*y/b) / (D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*(R^2*i^2*j^2));
        mkxy_0b(i,j) = vkxy;        mwxy_0b(i,j)=vwxy; % Insertar valores en la matriz
        mkxy_0b;
    end
    vkxy;
end
mkxy_0b;   mkxy_0b=sum(mkxy_0b,'all');      kxy_0b=mkxy_0b*(-32*R*a^2*p0)/(pi()^4)
mwxy_0b;   mwxy_0b=sum(mwxy_0b,'all');       wxy_0b=mwxy_0b*(16*p0*a^4)/(pi()^6)


for i = 1 : m % Coordenadas (x,y) = (a , b)
    vkxy=0;     vwxy=0;
    x=a;    y=b;
    for j = 1 : n
        vwxy = sin(i*pi()*x/a)*sin(j*pi()*y/b) / (i*j*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
        vkxy = cos(i*pi()*x/a)*cos(j*pi()*y/b) / (D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*(R^2*i^2*j^2));
        mkxy_ab(i,j) = vkxy;        mwxy_ab(i,j)=vwxy;% Insertar valores en la matriz
        mkxy_ab;
    end
    vkxy;
end
mkxy_ab;   mkxy_ab=sum(mkxy_ab,'all');      kxy_ab=mkxy_ab*(-32*R*a^2*p0)/(pi()^4)
mwxy_ab;   mwxy_ab=sum(mwxy_ab,'all');       wxy_ab=mwxy_ab*(16*p0*a^4)/(pi()^6)

disp('Donde desea saber w, kx, ky, kxy? ');% ESTE CALCULA PARA UN DETALLADO
x = input('X = ');
y = input('Y = ');
for i = 1: 2: m %Coordenadas (x,y) = (X,Y)
    i;    vwxy=0;    vkx = 0;  vky=0;
    for j = 1: 2 : n
        j;
        vwxy = sin(i*pi()*x/a)*sin(j*pi()*y/b) / (i*j*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
        vkx = i*sin(i*pi()*x/a)*sin(j*pi()*y/b) / (j*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
        vky = j*sin(i*pi()*x/a)*sin(j*pi()*y/b) / (i*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
        vkxy = cos(i*pi()*x/a)*cos(j*pi()*y/b) / (D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*(R^2*i^2*j^2));
        mwxy_c(i,j)=vwxy;     mkx(i,j) = vkx;       mky(i,j) = vky; mkxy_ab(i,j) = vkxy; % Insertar valores en la matriz
        mwxy_c;       mkx;        mky;
    end
    vwxy;   vkx;    vky;
end
mkxy_ab;   mkxy_ab=sum(mkxy_ab,'all');      kxy_D=mkxy_ab*(-32*R*a^2*p0)/(pi()^4)
mwxy_c;   mwxy_c=sum(mwxy_c,'all');      wxy_D=mwxy_c*(16*p0*a^4)/(pi()^6)
mkx;    mkx=sum(mkx,'all');              kxD=mkx*(16*p0*a^2)/(pi()^4)
mky;    mky=sum(mky,'all');              kyD=mky*(16*p0*a^2*R^2)/(pi()^4)


% ~~~   DEFORMACIONES   (Solo Curvaturas en este caso, por la naturaleza de la carga)
disp(' ');

disp('Coordenadas de análisis rápido: ');disp(' ');
disp('Caso 1: (a/2, b/2)');        
disp('Caso 2: (0, 0)');             
disp('Caso 3: (a, 0)');
disp('Caso 4: (0, b)');
disp('Caso 5: (a, b)');
disp(' ');
opcion = input('Seleccione el caso (n): ');

switch opcion % Selecciona las coordenadas de análisis rápido
    case 1
        msn_coordenadas='Usted eligio coordenadas (a/2, b/2)';
        disp (msn_coordenadas);
        ekxy=[0;0;0;    kx;   ky;  0];         % Coordenadas (a/2, b/2)
    case 2
        msn_coordenadas='Usted eligio coordenadas (0,0)';
        disp (msn_coordenadas);
        ekxy=[0;0;0;    0;    0;   kxy_00];    % Coordenadas (0, 0)
    case 3
        msn_coordenadas='Usted eligio coordenadas (a,0)';
        disp (msn_coordenadas);
        ekxy=[0;0;0;    0;    0;   kxy_a0];    % Coordenadas (a, 0)
    case 4
        msn_coordenadas='Usted eligio coordenadas (0,b)';
        disp (msn_coordenadas);
        ekxy=[0;0;0;    0;    0;   kxy_0b];    % Coordenadas (0, b)
    case 5
        msn_coordenadas='Usted eligio coordenadas (a,b)';
        disp (msn_coordenadas);
        ekxy=[0;0;0;    0;    0;   kxy_ab];    % Coordenadas (a, b)
end


ek=t;%        [mm]  Espesor de una capa/laminita
ec=ek*capas;%   [mm]  Espesor total del compuesto/estratificado
ec2cm =ec/10;   % [cm] espesor total del estratificado en cm
e_medio=ec/2;%  [mm]  Espesor medio, lugar de z=0 o plano medio
     
for k=1:capas               %       Cotas
    z(k)=-(e_medio-k*ek);%          Superiores
    z_1(k)=-(e_medio-(k-1)*ek);%    Inferiores
end

zk=z';              zk_1=z_1';

zks(ekxy,zk,zk_1,capas) % Realiza la multiplicación z{k}xy para superiores e inferiores
% Las variables de abajo son las que devuelve la función zks pero deben
% definirse globalmente para acceder a ellas fuera de la función
global kxsup; global kxinf; global kysup; global kyinf; global kxysup; global kxyinf

 % Recuerda solo tendras kurvaturas, no hay def en el plano ;) relajate
% -------------------------------------------------------------------
% 6.- Tensores de deformacion y esfuerzo globales y principales
% -------------------------------------------------------------------
% 7.- Aplicar teorías de falla  (°‿°)
% -------------------------------------------------------------------
% Asignar las deformaciones por cada capa /superiores /inferiores (for)
capas = length(orientaciones);
 disp('~~~~~~~~~~ TODO PARA CAPAS SUPERIORES ~~~~~~~~~~');
for i=1:capas % TODO PARA CAPAS SUPERIORES
      i;         
      % Defromaciones
      msncapa = ['Capa ',num2str(i),' s']; disp(msncapa);%Se descomenta
      exys = exyk_s(i,kxsup,kysup,kxysup);
            
      % Esfuerzos
      c = cosd(orientaciones(i));
      s = sind(orientaciones(i));
      
      % Transformación de la matriz de rigidez reducida osea Qxyk
      qactual = zeros(3,3); % Matriz vacia de 3x3 que se va llenando / se limpia

      qactual(1,1) = q11*c^4+q22*s^4+s*c*(s*c*(q21+q12+4*q33)+c^2*(2*q31+2*q13)+s^2*(2*q32+2*q23));
      qactual(1,2) = q12*c^4+q21*s^4+s*c*(s*c*(q11+q22-4*q33)+c^2*(2*q32-2*q13)+s^2*(2*q31-2*q23));
      qactual(1,3) = s*c*(s*c*(-2*q31+2*q32)+c^2*(-q11+q12)+s^2*(-q21+q22))+(c^2-s^2)*(q13*c^2+q23*s^2+2*q33*s*c);
      qactual(2,1) = q21*c^4+q12*s^4+s*c*(s*c*(q11+q22-4*q33)+c^2*(-2*q31+2*q23)+s^2*(-2*q32+2*q13));
      qactual(2,2) = q22*c^4+q11*s^4+s*c*(s*c*(q21+q12+4*q33)+c^2*(-2*q32-2*q23)+s^2*(-2*q31-2*q13));
      qactual(2,3) = s*c*(s*c*(2*q31-2*q32)+c^2*(-q21+q22)+s^2*(-q11+q12))+(c^2-s^2)*(q13*s^2+q23*c^2-2*q33*s*c);
      qactual(3,1) = s*c*(s*c*(-2*q13+2*q23)+c^2*(-q11+q21)+s^2*(-q12+q22)+2*q33*(c^2-s^2))+(c^2-s^2)*(q31*c^2+q32*s^2);
      qactual(3,2) = s*c*(s*c*(2*q13-2*q23)+c^2*(-q12+q22)+s^2*(-q11+q21)-2*q33*(c^2-s^2))+(c^2-s^2)*(q31*s^2+q32*c^2);
      qactual(3,3) = s*c*(s*c*(q11-q21-q12+q22)+(c^2-s^2)*(-q31+q32-q13+q23))+q33*(c^2-s^2)^2;
            
      sxys = qactual*exys; %    Esto da el esfuerzo de la capa actual
      
% Conversión A ejes principales para esfuerzos y deformaciones D'/ T'
        
      % Deformaciones en ejes pricncipales
        el  = exys(1)*c^2 + exys(2)*s^2 - exys(3)*s*c; % Deformación en l
        et  = exys(1)*s^2 + exys(2)*c^2 + exys(3)*s*c; % Deformación en t
        glt = 2*exys(1)*s*c - 2*exys(2)*s*c + exys(3)*(c^2-s^2); % Deformación angular gamma_lt
                
      % Esfuerzos en ejes principales
        sl  = sxys(1)*c^2 + sxys(2)*s^2 - 2*sxys(3)*s*c; % Sigma en l
        st  = sxys(1)*s^2 + sxys(2)*c^2 + 2*sxys(3)*s*c; % Sigma en t
        tlt = sxys(1)*s*c - sxys(2)*s*c + sxys(3)*(c^2-s^2); % Tau lt
        slt = [sl;st;tlt];
      
      % Aplicar una teoría de falla
        TsaiHillWu(slmaxT,stmaxT,slmaxC,stmaxC,taumax,slt);
        disp(' ');
        
end
disp('~~~~~~~~~~ TODO PARA CAPAS INFERIORES ~~~~~~~~~~');
for i=1:capas % TODO PARA CAPAS SUPERIORES
      i;         
      % Defromaciones
      msncapa = ['Capa ',num2str(i),' i']; disp(msncapa);%Se descomenta
      exyi = exyk_i(i,kxinf,kyinf,kxyinf);
            
      % Esfuerzos
      c = cosd(orientaciones(i));
      s = sind(orientaciones(i));
      
      % Transformación de la matriz de rigidez reducida osea Qxyk
      qactual = zeros(3,3); % Matriz vacia de 3x3 que se va llenando / se limpia

      qactual(1,1) = q11*c^4+q22*s^4+s*c*(s*c*(q21+q12+4*q33)+c^2*(2*q31+2*q13)+s^2*(2*q32+2*q23));
      qactual(1,2) = q12*c^4+q21*s^4+s*c*(s*c*(q11+q22-4*q33)+c^2*(2*q32-2*q13)+s^2*(2*q31-2*q23));
      qactual(1,3) = s*c*(s*c*(-2*q31+2*q32)+c^2*(-q11+q12)+s^2*(-q21+q22))+(c^2-s^2)*(q13*c^2+q23*s^2+2*q33*s*c);
      qactual(2,1) = q21*c^4+q12*s^4+s*c*(s*c*(q11+q22-4*q33)+c^2*(-2*q31+2*q23)+s^2*(-2*q32+2*q13));
      qactual(2,2) = q22*c^4+q11*s^4+s*c*(s*c*(q21+q12+4*q33)+c^2*(-2*q32-2*q23)+s^2*(-2*q31-2*q13));
      qactual(2,3) = s*c*(s*c*(2*q31-2*q32)+c^2*(-q21+q22)+s^2*(-q11+q12))+(c^2-s^2)*(q13*s^2+q23*c^2-2*q33*s*c);
      qactual(3,1) = s*c*(s*c*(-2*q13+2*q23)+c^2*(-q11+q21)+s^2*(-q12+q22)+2*q33*(c^2-s^2))+(c^2-s^2)*(q31*c^2+q32*s^2);
      qactual(3,2) = s*c*(s*c*(2*q13-2*q23)+c^2*(-q12+q22)+s^2*(-q11+q21)-2*q33*(c^2-s^2))+(c^2-s^2)*(q31*s^2+q32*c^2);
      qactual(3,3) = s*c*(s*c*(q11-q21-q12+q22)+(c^2-s^2)*(-q31+q32-q13+q23))+q33*(c^2-s^2)^2;
            
      sxyi = qactual*exyi; %    Esto da el esfuerzo de la capa actual
      
% Conversión A ejes principales para esfuerzos y deformaciones D'/ T'
        
      % Deformaciones en ejes pricncipales
        el  = exyi(1)*c^2 + exyi(2)*s^2 - exyi(3)*s*c; % Deformación en l
        et  = exyi(1)*s^2 + exyi(2)*c^2 + exyi(3)*s*c; % Deformación en t
        glt = 2*exyi(1)*s*c - 2*exyi(2)*s*c + exyi(3)*(c^2-s^2); % Deformación angular gamma_lt
        
      % Esfuerzos en ejes principales
        sl  = sxyi(1)*c^2 + sxyi(2)*s^2 - 2*sxyi(3)*s*c; % Sigma en l
        st  = sxyi(1)*s^2 + sxyi(2)*c^2 + 2*sxyi(3)*s*c; % Sigma en t
        tlt = sxyi(1)*s*c - sxyi(2)*s*c + sxyi(3)*(c^2-s^2); % Tau lt
        slt = [sl;st;tlt];
      
      % Aplicar una teoría de falla
        TsaiHillWu(slmaxT,stmaxT,slmaxC,stmaxC,taumax,slt);
        disp(' ');
        
end
disp('Haz terminado (☞ﾟ∀ﾟ)☞'); disp(' ');
disp(msn_coordenadas);
disp(' para factores de seguridad');
disp(' ')

disp('----------------------------------------------------------');
disp('CENTRO a/2,b/2');
centro = ['w = ',num2str(wxy_c),'| kx = ',num2str(kx),' | ky = ',num2str(ky)];
disp(centro); disp(' ');

disp('Borde 0,0');
borde00=['kxy = ',num2str(kxy_00),' | w = ',num2str(wxy_00)];
disp(borde00); disp(' ');

disp('Borde a,0');
bordea0=['kxy = ',num2str(kxy_a0),' | w = ',num2str(wxy_a0)];
disp(bordea0); disp(' ');

disp('Borde 0,b');
borde0b=['kxy = ',num2str(kxy_0b),' | w = ',num2str(wxy_0b)];
disp(borde0b); disp(' ');

disp('Borde a,b');
bordeab=['kxy = ',num2str(kxy_ab),' | w = ',num2str(wxy_ab)];
disp(bordeab);

disp(' ');msn_D = ['Coordenadas (',num2str(x),', ',num2str(y),')'];
disp(msn_D); valoresD2 = ['kx = ',num2str(kxD),' | ky = ',num2str(kyD)];
valoresD1 = ['w = ',num2str(wxy_D),' | kxy = ',num2str(kxy_D)];
disp(valoresD1);
disp(valoresD2);disp(' ');

%
i=0;j=0;
% _____________________ GRAFICANDO  ___________________________
mKx=[];  mKy=[];  mKxy=[];    mW=[];
mWxy=[];    mKxx=[];    mKyy=[];    mKxxyy=[];
for x=1:a% coordenada en x e y por mm
    for y=1:b
        
        for i=1:2:m
            for j=1:2:n
                vW = sin(i*pi()*x/a)*sin(j*pi()*y/b) / (i*j*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
                vKx = i*sin(i*pi()*x/a)*sin(j*pi()*y/b) / (j*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
                vKy = j*sin(i*pi()*x/a)*sin(j*pi()*y/b) / (i*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
                vKxy = cos(i*pi()*x/a)*cos(j*pi()*y/b) / (D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*(R^2*i^2*j^2));
                
                % (+) grafica pa riba, (-) grafia pa abajo
                mW(i,j) = -vW;      mKx(i,j) = -vKx;     mKy(i,j) = -vKy;
                mKxy(i,j) = -vKxy;
            end
        end
        mW=sum(mW,'all');   W = mW*(16*p0*a^4)/(pi()^6); % Obtengo un valor de w para coordenadas (x,y)
        mKx=sum(mKx,'all'); Kx = mKx*(16*p0*a^2)/(pi()^4);
        mKy=sum(mKy,'all'); Ky = mKy*(16*p0*a^2*R^2)/(pi()^4);
        mKxy=sum(mKxy,'all'); Kxy = mKxy*(-32*R*a^2*p0)/(pi()^4);
        
        mWxy(y,x)= W; % Agrega un valor a la matriz a plotear (xfila, ycolumna)
        mKxx(y,x)= Kx;
        mKyy(y,x)= Ky;
        mKxxyy(y,x)= Kxy;
    end
end

figure(1)
mesh(mWxy);
title('w(x,y)'); xlabel('Largo [mm]');       ylabel('Ancho [mm]');      zlabel('Espesor [mm]');

figure(2)
mesh(mKxx);
title('kx'); xlabel('Largo [mm]');       ylabel('Ancho [mm]');      zlabel('[1/mm]');

figure(3)
mesh(mKyy);
title('ky'); xlabel('Largo [mm]');       ylabel('Ancho [mm]');      zlabel('[1/mm]');

figure(4)
mesh(mKxxyy);
title('kxy'); xlabel('Largo [mm]');       ylabel('Ancho [mm]');      zlabel('[1/mm]');    
%}

% SACANDO LOS ESFUERZOS CORTANTES
% Para punto de análisis AQUI CAMBIA EL VALOR DE x y de y por eso hasta el
% final

disp('Coordenadas de analisis para cortantes:');
x = input('X = ');
y = input('Y = ');
coordenadas=['Para coordenadas: ','(',num2str(x),',' num2str(y),')'];
disp(coordenadas);
for i=1:2: m
    vQx_xy=0;
    vQy_xy=0;
    for j=1:2:n
        vQx_xy = (D11*i^2 + R^2*j^2*D12 + 2*R^2*j^2*D33)*(cos(i*pi()*x/a)*sin(j*pi()*y/b)) / (j*(D11*i^4 + (2*D12 + 4*D33)*(i^2*j^2*R^2) + D22*j^4*R^4));
        vQy_xy = (D21*i^2 + R^2*j^2*D22 + 2*R^2*j^2*D33)*(sin(i*pi()*x/a)*cos(j*pi()*y/b)) / (i*(D11*i^4 + (2*D12 + 4*D33)*(i^2*j^2*R^2) + D22*j^4*R^4));
        mQx_xy(i,j) = vQx_xy;
        mQy_xy(i,j) = vQy_xy;
    end
end
disp('Los Qx y Qy')
mQx_xy=sum(mQx_xy,'all');   Qx_xy = mQx_xy*(16*a*p0/(pi()^3));
mQy_xy=sum(mQy_xy,'all');   Qy_xy = mQy_xy*(16*a^2*p0/(b*pi()^3));
Qx_xy
Qy_xy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Encontrando las cargas en los apoyos :D

    % Igual ordenando por bordes 
  % Declarando matrices
  mFzQx=[]; mFzQy=[];
  
  disp('~~ CARGAS EN LOS APOYOS ~~');
  
  % Borde x=0 sobre y
  for i = 1:2: m
      vFzQx=0;
      for j = 1:2:n % n impar
          vFzQx = (D11*i^2 + R^2*j^2*D12 + 2*R^2*j^2*D33) / (j^2*(D11*i^4 + (2*D12 + 4*D33)*i^2*j^2*R^2 + D22*j^4*R^4));
          mFzQx(i,j) = vFzQx;
      end
  end
  mFzQx=sum(mFzQx,'all');   FzQx=mFzQx*(32*a*b*p0)/(pi()^4)

  % Borde y=0 sobre x
  for i = 1:2:m % m impar
      vFzQy=0;
      for j=1:2:n
          vFzQy = (D21*i^2 + D22*R^2*j^2 + 2*D33*i^2) / (i^2*(D11*i^4 + (2*D12 + 4*D33)*i^2*j^2*R^2 + D22*j^4*R^4));
          mFzQy(i,j) = vFzQy;
      end
  end
  mFzQy=sum(mFzQy,'all');   FzQy=mFzQy*(32*a^3*p0)/(b*pi()^4)

%   Calcularndo las CARGAS EN LOS APOYOS 

Fp0 = p0*a*b;
%   Carga en los bordes por cortantes
FzQ =  2*(FzQx + FzQy)

%   Fz en las esquinas
Fz_esquinas = Fp0 - FzQ

%   Fuerza en una esquina
Fz_esquina = Fz_esquinas/4

    
% ¡¡¡¡¡¡¡¡(☞ ಠ_ಠ)☞ ALTO APARTADO DE FUNCIONES  ¡¡¡¡¡¡¡¡¡¡
% Todo lo escrito debajo son funciones y deben ir al final del codigo
%..................................................................................
%        |||  TENSORES DE DEFORMACIONES GLOBALES  |||
    % Función 1.- Realiza la multiplicación z{k}xy para cada capa
    function [kxsup, kxinf, kysup, kyinf, kxysup , kxyinf] = zks(ekxy,zk,zk_1,capas)
    
         for p=1:capas
            %disp('Haz llegado al for')
            kx_sup(p) = ekxy(4)*zk(p);          kx_inf(p) = ekxy(4)*zk_1(p);
            ky_sup(p) = ekxy(5)*zk(p);          ky_inf(p) = ekxy(5)*zk_1(p);
            kxy_sup(p) = ekxy(6)*zk(p);         kxy_inf(p) = ekxy(6)*zk_1(p);
         end      
      %Reacomodo de las matrices
      global kxsup; global kxinf; global kysup; global kyinf; global kxysup;global kxyinf;
      kxsup=kx_sup';            kysup=ky_sup';          kxysup=kxy_sup';
      kxinf=kx_inf';            kyinf=ky_inf';          kxyinf=kxy_inf';  
    end
%...........................................................................    
%...........................................................................
 % Saca los tensores superiores de una capa dada
        %recordar que las deformaciones en membrana son 0
        %por eso {e}xy = {e⁰}+z{k};  es en este caso es {e}xy = z{k}
    function [exyk_s] = exyk_s(capa,kxs,kys,kxys) % Para las superiores
        exyk_s=[kxs(capa);kys(capa);kxys(capa)];
    end
    function [exyk_i] = exyk_i(capa,kxi,kyi,kxyi) % Para las inferiores
        exyk_i=[kxi(capa);kyi(capa);kxyi(capa)];
    end
%        |||  TENSORES DE ESFUERZOS GLOBALES  |||
    %los esfuersos son {s}xy = [Q]xy {e}xy
function [sxyk]=sxyk(qx,exyk)
    sxyk=qx*exyk;
end
%..........................................................................
% FACTORES DE RESISTENCIA por TSAI-HILL Y TSAI-WU

function [hill, Wu1, Wu2] = TsaiHillWu(slT,stT,slC,stC,taumax,slt)
    %Entradas
       
    sl = slt(1);
    st = slt(2);
    tau= slt(3);
    %Borrando calculos
    hill=0;
    s1T=0;
    s2T=0;
    s3=0;
    RTT=0;
    
    s1C=0;
    s2C=0;
    RCC=0;
    
    RTC=0;
    RCT=0;
    
    
    %~~~~~~~~~~ Por Tsai-Hill ~~~~~~~~~

disp('Tsai-Hill')
% Aqui va el discriminador de casos

if sl>0 && st>0 % Condición de TENSIÓN-TENSIÓN
    s1T=(sl^2-sl*st)/(slT^2);
    s2T=(st^2)/(stT^2);
    s3=(tau^2)/(taumax^2);
    %salida
    RTT=sqrt(1/(s1T+s2T+s3));
    hill=['RTT = ',num2str(RTT)];
    disp(hill)%hillTT
    
elseif sl<0 && st<0 %   COMPRESIÓN-COMPRESIÓN
    s1C=(sl^2-sl*st)/(slC^2);
    s2C=(st^2)/(stC^2);
    s3=(tau^2)/(taumax^2);
    %salida
    RCC=sqrt(1/(s1C+s2C+s3));
    hill=['RCC = ',num2str(RCC)];
    disp(hill)%hillCC
    
elseif sl>0 && st<0 %    TENSIÓN-COMPRESIÓN
    s1T=(sl^2-sl*st)/(slT^2);
    s2C=(st^2)/(stC^2);
    s3=(tau^2)/(taumax^2);
    %salida
    RTC=sqrt(1/(s1T+s2C+s3));
    hill=['RTC = ',num2str(RTC)];
    disp(hill)%hillTC
    
elseif sl<0 && st>0 %    COMPRESIÓN-TENSIÓN
    s1C=(sl^2-sl*st)/(slC^2);
    s2T=(st^2)/(stT^2);
    s3=(tau^2)/(taumax^2);
    %salida
    RCT=sqrt(1/(s1C+s2T+s3));
    hill=['RCT = ',num2str(RCT)];
    disp(hill)%hillCT
end

%   Borrando datos
f1=0;
f2=0;
f11=0;
f22=0;
f66=0;
f12=0;
a=0; b=0;
R1=0;   R2=0;
%~~~~~~~~~~ Por Tsai-Wu ~~~~~~~~~
%   Calculos
f1=1/slT-1/slC;
f2=1/stT-1/stC;
f11=1/(slT*slC);
f22=1/(stT*stC);
f66=1/(taumax^2);
f12=-sqrt(f11*f22)/2;

a=f11*sl^2+f22*st^2+f66*tau^2+2*f12*sl*st;
b=f1*sl+f2*st;
c=-1;
%   Salidas
disp('Tsai-Wu')
R1=(-b+sqrt(b^2-4*a*c))/(2*a);
R2=(-b-sqrt(b^2-4*a*c))/(2*a);
Wu1=['R1 = ',num2str(R1)];
Wu2=['R2 = ',num2str(R2)];
disp(Wu1)
disp(Wu2)   
end