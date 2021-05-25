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
## @deftypefn {} {@var{retval} =} zks (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Felipe Carmelo <chetos@chetos-Inspiron-5558>
## Created: 2021-05-02

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
