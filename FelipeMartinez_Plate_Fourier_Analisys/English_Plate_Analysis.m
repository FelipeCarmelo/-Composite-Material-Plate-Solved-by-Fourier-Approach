%{
´◔‿ゝ◔`)━☞ANALISYS OF A PLATE SIMPLY SUPPORTED IN ITS BORDERS AND LOADED IN THE CENTRE

************************* DESCRIPTION *****************************

Theres a plate of dimension "a" by "b" and loaded in its centre, the
domain of the load is a square surface of "c" dimensions, the load is in
Newtons and the dimensions in millimeters.
 
AND YOU CAN DEFINE

1) Number of layers
2) Orientations
3) The Material
4) The load

YOU OBTAIN 

1) 3D Graphics for axial & angular deformations and kurvatures
2) Tables of deformations/loads/shears at the poit of the plate
   domain that you define and in its borders by default
3) Two failure theories applied for each layer

    SO THIS IS LIKE A MODIFY AND SEE WHAT HAPPEN
    
                    HAVE FUN !!!
*******************************************************************
 %}
clear;close all; clear all; clc;

% Orientations
orientations = [60 -60 -60 60 60 -60 -60 60 60 -60 -60 60];

% Plate Dimensions
long = 800; % mm "a"
width = 600;  % mm "b"

% Number of layers
layers = length(orientations);

%__________________________________________________________________________
%         MECHANICAL PROPERTIES OF MATERIAL (l,t) 

% Glass VE BY DEFAULT
         El = 31161   ;% MPA
         Et = 7452   ; % MPA
        Glt = 2738   ; % MPA
        vlt = 0.3397   ;
        vtl = vlt*(Et/El);  
        
        slmaxT = 500   ;%   MPA
        slmaxC = 500   ;%   MPA
        stmaxT = 50   ; %   MPA
        stmaxC = 50   ; %   MPA
        taumax = 50   ; %     MPA
        
        % Thickness
        
        t=0.2; % [ mm ]

% 1.- Finding [S]lt  &   [Q]lt
s11=1/El;   s12=-vtl/Et;    s21=-vlt/El;    s22=1/Et;   s66=1/Glt;

Slt=[s11    s12     0;
     s21    s22     0;
     0      0       s66];
 
Qlt = inv(Slt);

q11 = Qlt(1,1);     q12 = Qlt(1,2);     q13 = Qlt(1,3);
q21 = Qlt(2,1);     q22 = Qlt(2,2);     q23 = Qlt(2,3);
q31 = Qlt(3,1);     q32 = Qlt(3,2);     q33 = Qlt(3,3);

% -------------------------------------------------------------------
%                 ABD & abd matrices 

% Some Stuff that let the matrices work
% TAKE CARE OF TOUCHING THESE LINES
N = layers; h = zeros (N+1); h(N+1) = N*t/2;   qbarK = zeros(3*N,3);       

for k = 1 : N
    c = cosd(orientations(k));
    s = sind(orientations(k));
    
    % Stiffness matrix transformation
    qactual = zeros(3,3); % Container matrix/it takes values of a matrix by the layer orientation and then the next layer and so on/ it overwrites itself by each cycle
    
    qactual(1,1) = q11*c^4+q22*s^4+s*c*(s*c*(q21+q12+4*q33)+c^2*(2*q31+2*q13)+s^2*(2*q32+2*q23));
    qactual(1,2) = q12*c^4+q21*s^4+s*c*(s*c*(q11+q22-4*q33)+c^2*(2*q32-2*q13)+s^2*(2*q31-2*q23));
    qactual(1,3) = s*c*(s*c*(-2*q31+2*q32)+c^2*(-q11+q12)+s^2*(-q21+q22))+(c^2-s^2)*(q13*c^2+q23*s^2+2*q33*s*c);
    qactual(2,1) = q21*c^4+q12*s^4+s*c*(s*c*(q11+q22-4*q33)+c^2*(-2*q31+2*q23)+s^2*(-2*q32+2*q13));
    qactual(2,2) = q22*c^4+q11*s^4+s*c*(s*c*(q21+q12+4*q33)+c^2*(-2*q32-2*q23)+s^2*(-2*q31-2*q13));
    qactual(2,3) = s*c*(s*c*(2*q31-2*q32)+c^2*(-q21+q22)+s^2*(-q11+q12))+(c^2-s^2)*(q13*s^2+q23*c^2-2*q33*s*c);
    qactual(3,1) = s*c*(s*c*(-2*q13+2*q23)+c^2*(-q11+q21)+s^2*(-q12+q22)+2*q33*(c^2-s^2))+(c^2-s^2)*(q31*c^2+q32*s^2);
    qactual(3,2) = s*c*(s*c*(2*q13-2*q23)+c^2*(-q12+q22)+s^2*(-q11+q21)-2*q33*(c^2-s^2))+(c^2-s^2)*(q31*s^2+q32*c^2);
    qactual(3,3) = s*c*(s*c*(q11-q21-q12+q22)+(c^2-s^2)*(-q31+q32-q13+q23))+q33*(c^2-s^2)^2;
   
    qbarK([3*k-2:3*k],[1:3])=qactual([1:3],[1:3]); % Matrix that keeps all the qactual's cycling matrices
    h(k)=(k-N/2-1)*t;   % Values of the positions of each layer in the overall composite thickness

end
A = zeros(3);   B = zeros(3);   D = zeros(3);


for i=1:3 % Cycle that make the products sum of(Qij)(zk-zk1)² | (Qij)(zk-zk1)³
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

    if rem(length(orientations),2)==0
    A= A.*[1 1 0 ; 1 1 0 ; 0 0 1];
    B= B.*[0 0 0 ; 0 0 0 ; 0 0 0];
    D= D.*[1 1 1 ; 1 1 1 ; 1 1 1];
    end
    A
    B
    D
% ABD and abd matrices
ABD = [A,B;B,D];                    abd = inv(ABD);

% Isolated Dij elements
D11 = D(1,1);   D12 = D(1,2);   D13 = D(1,3);
D21 = D(2,1);   D22 = D(2,2);   D23 = D(2,3);
D31 = D(3,1);   D32 = D(3,2);   D33 = D(3,3);

% -------------------------------------------------------------------
% Calculations of kurvatures, kx, ky, kxy and w(x,y) deformation

    % LOAD DEFINITION
gravity = 9.81; % m/s2
Mass = 0.2; % kg
load_surface = 200*200; % mm2 "c"x"c"
p0=Mass*gravity/load_surface;% MPa

 % ~~~~~~~~~~~ FINDING w(x,y), kx, ky, kxy ~~~~~~~~~~~~
        % Common Coeficients
    a=long;    b=width;    R=a/b;  
        
        % Number of iterations
    mn=input('Define the iterations number: ');
    m=2*mn;   % m = rows     (i)
    n=m;    % n = colums  (j)
    
                % Initial Matrices (empties)
    mwxy_c=[];      mwxy_00=[];     mwxy_a0=[];     mwxy_0b=[];     mwxy_ab=[]; 
    mkx=[];     mky=[];    
    mkxy_00=[];     mkxy_a0=[];     mkxy_0b=[];     mkxy_ab=[];     
    
    % *******************************************************
    % Next the calculations shown in the images are performed
    % *******************************************************
    
    %   For kx, ky y w(x,y) in the centre
for i = 1: 2: m % Coordinates(x,y) = (a/2 ,b/2) with step of  2 to take odd numbers
    i;    vwxy=0;    vkx = 0;  vky=0;
    x=a/2;    y=b/2;
    for j = 1: 2 : n
        j;
        vwxy = sin(i*pi()*x/a)*sin(j*pi()*y/b) / (i*j*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
        vkx = i*sin(i*pi()*x/a)*sin(j*pi()*y/b) / (j*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
        vky = j*sin(i*pi()*x/a)*sin(j*pi()*y/b) / (i*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
        mwxy_c(i,j)=vwxy;     mkx(i,j) = vkx;       mky(i,j) = vky;
        mwxy_c;       mkx;        mky;
    end
    vwxy;   vkx;    vky;
end
mwxy_c;   mwxy_c=sum(mwxy_c,'all');      wxy_c=mwxy_c*(16*p0*a^4)/(pi()^6)
mkx;    mkx=sum(mkx,'all');              kx=mkx*(16*p0*a^2)/(pi()^4)
mky;    mky=sum(mky,'all');              ky=mky*(16*p0*a^2*R^2)/(pi()^4)

%   For w(x,y) and kxy in the borders
for i = 1 : m % Coordinates (x,y) = (0 , 0)
    vkxy=0; vwxy=0;
    x=0;    y=0;
    for j = 1 : n
        vwxy = sin(i*pi()*x/a)*sin(j*pi()*y/b) / (i*j*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
        vkxy = cos(i*pi()*x/a)*cos(j*pi()*y/b) / (D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*(R^2*i^2*j^2));
        mkxy_00(i,j) = vkxy;        mwxy_00(i,j)=vwxy; 
        mkxy_00;
    end
    vkxy;
end
mkxy_00;   mkxy_00=sum(mkxy_00,'all');      kxy_00=mkxy_00*(-32*R*a^2*p0)/(pi()^4)
mwxy_00;   mwxy_00=sum(mwxy_00,'all');      wxy_00=mwxy_00*(16*p0*a^4)/(pi()^6)


for i = 1 : m % Coordinates(x,y) = (a , 0)
    vkxy=0;     vwxy=0;
    x=a;    y=0;
    for j = 1 : n
        vwxy = sin(i*pi()*x/a)*sin(j*pi()*y/b) / (i*j*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
        vkxy = cos(i*pi()*x/a)*cos(j*pi()*y/b) / (D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*(R^2*i^2*j^2));
        mkxy_a0(i,j) = vkxy;    mwxy_a0(i,j)=vwxy;
        mkxy_a0;
    end
    vkxy;
end
mkxy_a0;   mkxy_a0=sum(mkxy_a0,'all');      kxy_a0=mkxy_a0*(-32*R*a^2*p0)/(pi()^4)
mwxy_a0;   mwxy_a0=sum(mwxy_a0,'all');       wxy_a0=mwxy_a0*(16*p0*a^4)/(pi()^6)


for i = 1 : m % Coordinates (x,y) = (0 , b)
    vkxy=0;     vwxy=0;
    x=0;    y=b;
    for j = 1 : n
        vwxy = sin(i*pi()*x/a)*sin(j*pi()*y/b) / (i*j*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
        vkxy = cos(i*pi()*x/a)*cos(j*pi()*y/b) / (D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*(R^2*i^2*j^2));
        mkxy_0b(i,j) = vkxy;        mwxy_0b(i,j)=vwxy;
        mkxy_0b;
    end
    vkxy;
end
mkxy_0b;   mkxy_0b=sum(mkxy_0b,'all');      kxy_0b=mkxy_0b*(-32*R*a^2*p0)/(pi()^4)
mwxy_0b;   mwxy_0b=sum(mwxy_0b,'all');       wxy_0b=mwxy_0b*(16*p0*a^4)/(pi()^6)


for i = 1 : m % Coordinates (x,y) = (a , b)
    vkxy=0;     vwxy=0;
    x=a;    y=b;
    for j = 1 : n
        vwxy = sin(i*pi()*x/a)*sin(j*pi()*y/b) / (i*j*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
        vkxy = cos(i*pi()*x/a)*cos(j*pi()*y/b) / (D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*(R^2*i^2*j^2));
        mkxy_ab(i,j) = vkxy;        mwxy_ab(i,j)=vwxy;
        mkxy_ab;
    end
    vkxy;
end
mkxy_ab;   mkxy_ab=sum(mkxy_ab,'all');      kxy_ab=mkxy_ab*(-32*R*a^2*p0)/(pi()^4)
mwxy_ab;   mwxy_ab=sum(mwxy_ab,'all');       wxy_ab=mwxy_ab*(16*p0*a^4)/(pi()^6)

disp('In wich coordinates would you like to knwo w, kx, ky, kxy? ');% ESTE CALCULA PARA UN DETALLADO
x = input('X = ');
y = input('Y = ');
for i = 1: 2: m % ASKED Coordinates (x,y) = (X,Y)
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


% ~~~   Strains   (Only Kurvatures in this case by the nature of the load)
disp(' ');

disp('Quick analisys coordinates: ');disp(' ');
disp('Case 1: (a/2, b/2)');        
disp('Case 2: (0, 0)');             
disp('Case 3: (a, 0)');
disp('Case 4: (0, b)');
disp('Case 5: (a, b)');
disp(' ');
opcion = input('Choose the case (n): ');

switch opcion 
    case 1
        msn_coordenadas='You choosed (a/2, b/2)';
        disp (msn_coordenadas);
        ekxy=[0;0;0;    kx;   ky;  0];        
    case 2
        msn_coordenadas='You choosed (0,0)';
        disp (msn_coordenadas);
        ekxy=[0;0;0;    0;    0;   kxy_00];  
    case 3
        msn_coordenadas='You choosed (a,0)';
        disp (msn_coordenadas);
        ekxy=[0;0;0;    0;    0;   kxy_a0];    
    case 4
        msn_coordenadas='You choosed (0,b)';
        disp (msn_coordenadas);
        ekxy=[0;0;0;    0;    0;   kxy_0b];   
    case 5
        msn_coordenadas='You choosed (a,b)';
        disp (msn_coordenadas);
        ekxy=[0;0;0;    0;    0;   kxy_ab];    
end


ek=t;%        [mm]  Thickness of a layer
ec=ek*layers;%   [mm]  Thickness of the overall plate (all layers)
ec2cm =ec/10;   % [cm] Thickness
e_medio=ec/2;%  [mm]  Half plate thickness, place where z=0 or mid plane
     
for k=1:layers               %      Positions of :
    z(k)=-(e_medio-k*ek);%          Superior faces
    z_1(k)=-(e_medio-(k-1)*ek);%    Lower faces
end

zk=z';              zk_1=z_1';

zks(ekxy,zk,zk_1,layers) % A function that obtains properties about kurvatures of the laminate
global kxsup; global kxinf; global kysup; global kyinf; global kxysup; global kxyinf

% Theres no plain strain due load nature

% Then we assign the strains by each layer face (superior or lower)
%layers = length(orientations);
 disp('~~~~~~~~~~ For superior faces ~~~~~~~~~~');
for i=1:layers
      i;         
      % Strains
      msncapa = ['Layer ',num2str(i),' s']; disp(msncapa);
      exys = exyk_s(i,kxsup,kysup,kxysup);
            
      % Stresses
      c = cosd(orientations(i));
      s = sind(orientations(i));
      
      % Stiffness matrix transformation Qxyk to GLOBAL AXES
      qactual = zeros(3,3); 

      qactual(1,1) = q11*c^4+q22*s^4+s*c*(s*c*(q21+q12+4*q33)+c^2*(2*q31+2*q13)+s^2*(2*q32+2*q23));
      qactual(1,2) = q12*c^4+q21*s^4+s*c*(s*c*(q11+q22-4*q33)+c^2*(2*q32-2*q13)+s^2*(2*q31-2*q23));
      qactual(1,3) = s*c*(s*c*(-2*q31+2*q32)+c^2*(-q11+q12)+s^2*(-q21+q22))+(c^2-s^2)*(q13*c^2+q23*s^2+2*q33*s*c);
      qactual(2,1) = q21*c^4+q12*s^4+s*c*(s*c*(q11+q22-4*q33)+c^2*(-2*q31+2*q23)+s^2*(-2*q32+2*q13));
      qactual(2,2) = q22*c^4+q11*s^4+s*c*(s*c*(q21+q12+4*q33)+c^2*(-2*q32-2*q23)+s^2*(-2*q31-2*q13));
      qactual(2,3) = s*c*(s*c*(2*q31-2*q32)+c^2*(-q21+q22)+s^2*(-q11+q12))+(c^2-s^2)*(q13*s^2+q23*c^2-2*q33*s*c);
      qactual(3,1) = s*c*(s*c*(-2*q13+2*q23)+c^2*(-q11+q21)+s^2*(-q12+q22)+2*q33*(c^2-s^2))+(c^2-s^2)*(q31*c^2+q32*s^2);
      qactual(3,2) = s*c*(s*c*(2*q13-2*q23)+c^2*(-q12+q22)+s^2*(-q11+q21)-2*q33*(c^2-s^2))+(c^2-s^2)*(q31*s^2+q32*c^2);
      qactual(3,3) = s*c*(s*c*(q11-q21-q12+q22)+(c^2-s^2)*(-q31+q32-q13+q23))+q33*(c^2-s^2)^2;
            
      sxys = qactual*exys; %    Kepts the current layer stress
      
% Strain and stresses transformation D' and T' to principal axes
        
      % Strains in principal axes
        el  = exys(1)*c^2 + exys(2)*s^2 - exys(3)*s*c; % longitudinal strain l
        et  = exys(1)*s^2 + exys(2)*c^2 + exys(3)*s*c; % transversal strain t
        glt = 2*exys(1)*s*c - 2*exys(2)*s*c + exys(3)*(c^2-s^2); % angular strain gamma_lt
        
      % Stresses in principal axes
        sl  = sxys(1)*c^2 + sxys(2)*s^2 - 2*sxys(3)*s*c; % Sigma l
        st  = sxys(1)*s^2 + sxys(2)*c^2 + 2*sxys(3)*s*c; % Sigma t
        tlt = sxys(1)*s*c - sxys(2)*s*c + sxys(3)*(c^2-s^2); % Tau lt
        slt = [sl;st;tlt];
      
      % Applied failure theory
        TsaiHillWu(slmaxT,stmaxT,slmaxC,stmaxC,taumax,slt);
        disp(' ');
        
end
disp('~~~~~~~~~~ For lower faces ~~~~~~~~~~');
for i=1:layers 
      i;         
      % Strains
      msncapa = ['Layer ',num2str(i),' i']; disp(msncapa);
      exyi = exyk_i(i,kxinf,kyinf,kxyinf);
            
      % Stresses
      c = cosd(orientations(i));
      s = sind(orientations(i));
      
      % Stiffness matrix transformation Qxyk
      qactual = zeros(3,3);

      qactual(1,1) = q11*c^4+q22*s^4+s*c*(s*c*(q21+q12+4*q33)+c^2*(2*q31+2*q13)+s^2*(2*q32+2*q23));
      qactual(1,2) = q12*c^4+q21*s^4+s*c*(s*c*(q11+q22-4*q33)+c^2*(2*q32-2*q13)+s^2*(2*q31-2*q23));
      qactual(1,3) = s*c*(s*c*(-2*q31+2*q32)+c^2*(-q11+q12)+s^2*(-q21+q22))+(c^2-s^2)*(q13*c^2+q23*s^2+2*q33*s*c);
      qactual(2,1) = q21*c^4+q12*s^4+s*c*(s*c*(q11+q22-4*q33)+c^2*(-2*q31+2*q23)+s^2*(-2*q32+2*q13));
      qactual(2,2) = q22*c^4+q11*s^4+s*c*(s*c*(q21+q12+4*q33)+c^2*(-2*q32-2*q23)+s^2*(-2*q31-2*q13));
      qactual(2,3) = s*c*(s*c*(2*q31-2*q32)+c^2*(-q21+q22)+s^2*(-q11+q12))+(c^2-s^2)*(q13*s^2+q23*c^2-2*q33*s*c);
      qactual(3,1) = s*c*(s*c*(-2*q13+2*q23)+c^2*(-q11+q21)+s^2*(-q12+q22)+2*q33*(c^2-s^2))+(c^2-s^2)*(q31*c^2+q32*s^2);
      qactual(3,2) = s*c*(s*c*(2*q13-2*q23)+c^2*(-q12+q22)+s^2*(-q11+q21)-2*q33*(c^2-s^2))+(c^2-s^2)*(q31*s^2+q32*c^2);
      qactual(3,3) = s*c*(s*c*(q11-q21-q12+q22)+(c^2-s^2)*(-q31+q32-q13+q23))+q33*(c^2-s^2)^2;
            
      sxyi = qactual*exyi;
      
% Strain and stresses transformation D' and T' to principal axes
        
       % Strains in principal axes
        el  = exyi(1)*c^2 + exyi(2)*s^2 - exyi(3)*s*c; % Deformación en l
        et  = exyi(1)*s^2 + exyi(2)*c^2 + exyi(3)*s*c; % Deformación en t
        glt = 2*exyi(1)*s*c - 2*exyi(2)*s*c + exyi(3)*(c^2-s^2); % Deformación angular gamma_lt
        
       % Stresses in principal axes
        sl  = sxyi(1)*c^2 + sxyi(2)*s^2 - 2*sxyi(3)*s*c; % Sigma en l
        st  = sxyi(1)*s^2 + sxyi(2)*c^2 + 2*sxyi(3)*s*c; % Sigma en t
        tlt = sxyi(1)*s*c - sxyi(2)*s*c + sxyi(3)*(c^2-s^2); % Tau lt
        slt = [sl;st;tlt];
      
       % Applied failure theory
        TsaiHillWu(slmaxT,stmaxT,slmaxC,stmaxC,taumax,slt);
        disp(' ');
        
end
disp('DONE (☞ﾟ∀ﾟ)☞'); disp(' ');
disp(msn_coordenadas);
disp(' Failure Theory');
%disp (msn_finish);
disp(' ')

disp('----------------------------------------------------------');
disp('CENTRE a/2,b/2');
centro = ['w = ',num2str(wxy_c),'| kx = ',num2str(kx),' | ky = ',num2str(ky)];
disp(centro); disp(' ');

disp('Border 0,0');
borde00=['kxy = ',num2str(kxy_00),' | w = ',num2str(wxy_00)];
disp(borde00); disp(' ');

disp('Border a,0');
bordea0=['kxy = ',num2str(kxy_a0),' | w = ',num2str(wxy_a0)];
disp(bordea0); disp(' ');

disp('Border 0,b');
borde0b=['kxy = ',num2str(kxy_0b),' | w = ',num2str(wxy_0b)];
disp(borde0b); disp(' ');

disp('Border a,b');
bordeab=['kxy = ',num2str(kxy_ab),' | w = ',num2str(wxy_ab)];
disp(bordeab);

disp(' ');msn_D = ['Coordinates (',num2str(x),', ',num2str(y),')'];
disp(msn_D); valoresD2 = ['kx = ',num2str(kxD),' | ky = ',num2str(kyD)];
valoresD1 = ['w = ',num2str(wxy_D),' | kxy = ',num2str(kxy_D)];
disp(valoresD1);
disp(valoresD2);disp(' ');

%
i=0;j=0;

% _______________________   G R A P H I C S  __________________________________
mKx=[];  mKy=[];  mKxy=[];    mW=[];
mWxy=[];    mKxx=[];    mKyy=[];    mKxxyy=[];
for x=1:a  % coordinates x and y by millimeter
    for y=1:b
        
        for i=1:2:m
            for j=1:2:n
                vW = sin(i*pi()*x/a)*sin(j*pi()*y/b) / (i*j*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
                vKx = i*sin(i*pi()*x/a)*sin(j*pi()*y/b) / (j*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
                vKy = j*sin(i*pi()*x/a)*sin(j*pi()*y/b) / (i*(D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*((R^2)*(i^2)*(j^2))));
                vKxy = cos(i*pi()*x/a)*cos(j*pi()*y/b) / (D11*(i^4)+D22*(R^4)*(j^4)+(2*D12+4*D33)*(R^2*i^2*j^2));
            
                mW(i,j) = -vW;      mKx(i,j) = -vKx;     mKy(i,j) = -vKy;
                mKxy(i,j) = -vKxy;
            end
        end
        mW=sum(mW,'all');   W = mW*(16*p0*a^4)/(pi()^6);
        mKx=sum(mKx,'all'); Kx = mKx*(16*p0*a^2)/(pi()^4);
        mKy=sum(mKy,'all'); Ky = mKy*(16*p0*a^2*R^2)/(pi()^4);
        mKxy=sum(mKxy,'all'); Kxy = mKxy*(-32*R*a^2*p0)/(pi()^4);
        
        % Adds a value for the matrix that plots (x=row, y=column)
        mWxy(y,x)= W; 
        mKxx(y,x)= Kx;
        mKyy(y,x)= Ky;
        mKxxyy(y,x)= Kxy;
    end
end

figure(1)
mesh(mWxy);
title('w(x,y)'); xlabel('Length [mm]');       ylabel('Width [mm]');      zlabel('Thickness [mm]');

figure(2)
mesh(mKxx);
title('kx'); xlabel('Length [mm]');       ylabel('Width [mm]');      zlabel('[1/mm]');

figure(3)
mesh(mKyy);
title('ky'); xlabel('Length [mm]');       ylabel('Width [mm]');      zlabel('[1/mm]');

figure(4)
mesh(mKxxyy);
title('kxy'); xlabel('Length [mm]');       ylabel('Width [mm]');      zlabel('[1/mm]');

% Finding shears stresses 

disp('Coordinates of analisys for shears:');
x = input('X = ');
y = input('Y = ');
coordinates=['For coordinates: ','(',num2str(x),',' num2str(y),')'];
disp(coordinates);
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
disp('The Qx and Qy')
mQx_xy=sum(mQx_xy,'all');   Qx_xy = mQx_xy*(16*a*p0/(pi()^3));
mQy_xy=sum(mQy_xy,'all');   Qy_xy = mQy_xy*(16*a^2*p0/(b*pi()^3));
Qx_xy
Qy_xy


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Reactions in the supports :D

  mFzQx=[]; mFzQy=[];
  
  disp('~~ LOADS IN THE SUPPORTS ~~');
  
  % Border x=0, along y
  for i = 1:2: m
      vFzQx=0;
      for j = 1:2:n % n odd
          vFzQx = (D11*i^2 + R^2*j^2*D12 + 2*R^2*j^2*D33) / (j^2*(D11*i^4 + (2*D12 + 4*D33)*i^2*j^2*R^2 + D22*j^4*R^4));
          mFzQx(i,j) = vFzQx;
      end
  end
  mFzQx=sum(mFzQx,'all');   FzQx=mFzQx*(32*a*b*p0)/(pi()^4)

  % Border y=0, along x
  for i = 1:2:m % m impar
      vFzQy=0;
      for j=1:2:n
          vFzQy = (D21*i^2 + D22*R^2*j^2 + 2*D33*i^2) / (i^2*(D11*i^4 + (2*D12 + 4*D33)*i^2*j^2*R^2 + D22*j^4*R^4));
          mFzQy(i,j) = vFzQy;
      end
  end
  mFzQy=sum(mFzQy,'all');   FzQy=mFzQy*(32*a^3*p0)/(b*pi()^4)

Fp0 = p0*a*b;
%   Shear loads in borders 
FzQ =  2*(FzQx + FzQy)

%   Frocers in the edges Fz
Fz_edges = Fp0 - FzQ

%   
Fz_edge = Fz_edges/4

    
% ¡¡¡¡¡¡¡¡(☞ ಠ_ಠ)☞ STOP, FUNCTIONS ANNEX  ¡¡¡¡¡¡¡¡¡¡
% All the next lines must be at the end of the code
%..................................................................................
%        |||  GLOBAL STRAINS TENSORS  |||
    % Performs products z{k}xy for each layer
    function [kxsup, kxinf, kysup, kyinf, kxysup , kxyinf] = zks(ekxy,zk,zk_1,layers)
    
         for p=1:layers
            kx_sup(p) = ekxy(4)*zk(p);          kx_inf(p) = ekxy(4)*zk_1(p);
            ky_sup(p) = ekxy(5)*zk(p);          ky_inf(p) = ekxy(5)*zk_1(p);
            kxy_sup(p) = ekxy(6)*zk(p);         kxy_inf(p) = ekxy(6)*zk_1(p);
         end      
      % Rearanging matrices
      global kxsup; global kxinf; global kysup; global kyinf; global kxysup;global kxyinf;
      kxsup=kx_sup';            kysup=ky_sup';          kxysup=kxy_sup';
      kxinf=kx_inf';            kyinf=ky_inf';          kxyinf=kxy_inf';  
    end 
%...........................................................................

 % Get tensors for uper and lower faces
        % Theres only {e}xy = {e⁰}+z{k};  in this case {e}xy = z{k}
    function [exyk_s] = exyk_s(capa,kxs,kys,kxys) % Upers faces
        exyk_s=[kxs(capa);kys(capa);kxys(capa)];
    end
    function [exyk_i] = exyk_i(capa,kxi,kyi,kxyi) % Lower faces
        exyk_i=[kxi(capa);kyi(capa);kxyi(capa)];
    end
%        |||  GLOBAL STRESSES TENSORS   |||
    % Strssses are {s}xy = [Q]xy {e}xy
function [sxyk]=sxyk(qx,exyk)
    sxyk=qx*exyk;
end
%..........................................................................

% FACTORES DE RESISTENCIA 

% BY TSAI-HILL and TSAI-WU CRITERIA

function [hill, Wu1, Wu2] = TsaiHillWu(slT,stT,slC,stC,taumax,slt)
    % Inputs
       
    sl = slt(1);
    st = slt(2);
    tau= slt(3);
    % Cleaning old variables
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
    
    
    %~~~~~~~~~~ Tsai-Hill Criterion ~~~~~~~~~

disp('Tsai-Hill')

% Programming different conditions
if sl>0 && st>0 % TENSION-TENSION
    s1T=(sl^2-sl*st)/(slT^2);
    s2T=(st^2)/(stT^2);
    s3=(tau^2)/(taumax^2);
    % output
    RTT=sqrt(1/(s1T+s2T+s3));
    hill=['RTT = ',num2str(RTT)];
    disp(hill)%hillTT
    
elseif sl<0 && st<0 %   Compress-compress
    s1C=(sl^2-sl*st)/(slC^2);
    s2C=(st^2)/(stC^2);
    s3=(tau^2)/(taumax^2);
    % out
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
    
elseif sl<0 && st>0 %    COMPRESS-TENSION
    s1C=(sl^2-sl*st)/(slC^2);
    s2T=(st^2)/(stT^2);
    s3=(tau^2)/(taumax^2);
    % Out
    RCT=sqrt(1/(s1C+s2T+s3));
    hill=['RCT = ',num2str(RCT)];
    disp(hill)%hillCT
end

%   Cleaning coeficients
f1=0;
f2=0;
f11=0;
f22=0;
f66=0;
f12=0;
a=0; b=0;
R1=0;   R2=0;
%~~~~~~~~~~ Tsai-Wu Criterion ~~~~~~~~~
%   Formulas
f1=1/slT-1/slC;
f2=1/stT-1/stC;
f11=1/(slT*slC);
f22=1/(stT*stC);
f66=1/(taumax^2);
f12=-sqrt(f11*f22)/2;

a=f11*sl^2+f22*st^2+f66*tau^2+2*f12*sl*st;
b=f1*sl+f2*st;
c=-1;
%   Outputs
disp('Tsai-Wu')
R1=(-b+sqrt(b^2-4*a*c))/(2*a);
R2=(-b-sqrt(b^2-4*a*c))/(2*a);
Wu1=['R1 = ',num2str(R1)];
Wu2=['R2 = ',num2str(R2)];
disp(Wu1)
disp(Wu2)   
end