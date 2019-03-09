 function [ CDF ] = Davies_inversion_continuous( CF, my_range )
%function [ CDF ] = Davies_inversion_continuous( CF, my_range )
%   Numerically inverts the Characteristic function CF (sym object) over the range
%   my_range. It outputs the continuous CDF. symvar(CF) must not be called x. 
%   INPUT:
%       - CF: sym object representing the expression of the characteristic
%       function. symvar(CF) must not be called x. Must have only one
%       unknown.
%       - my_range: array of points where to calculate the CDF of CF
%   OUTPUT:
%       - CDF: the value of the CDF at the points specified in my_range
%
%   Inversion
%   according to Davies, Robert B. "Numerical inversion of a characteristic function." Biometrika 60.2 (1973): 415-417.
%   
%   Copyright (C) 2017  Antonio Franco
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

CDF = zeros(1, size(my_range,2));

t = symvar(CF);

syms 'x' real;

Integrand = imag(CF*exp(-1i*t*x)/(2*pi*t));

for i =1:numel(my_range)
    xx = my_range(i);
    my_integrand = subs(Integrand, x, xx);
    try
        my_integral = double(vpaintegral(my_integrand, t, -inf, inf));
    catch
        XX = linspace(-500,500,1e4);
        fun = matlabFunction(my_integrand);
        my_integral = trapz(XX,fun(XX));
    end
    CDF(i) = 0.5 - my_integral;
end

end