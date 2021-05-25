## Copyright (C) 2021 Felipe Carmelo
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {} {@var{retval} =} TsaiHillWu (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Felipe Carmelo <chetos@chetos-Inspiron-5558>
## Created: 2021-05-02

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