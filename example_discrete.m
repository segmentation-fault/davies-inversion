%Example on how to use Davies_inversion_discrete.m

clc
clear

%Discrete uniform
syms 't'
a = 3;
b = 10;
n = b-a+1;
CF1 = (exp(1i*a*t) - exp(1i*(b+1)*t))/(n*(1-exp(1i*t)));
my_range = a-2:b+2;
CDF = Davies_inversion_discrete( CF1, my_range );
figure
subplot(2,2,1);
stairs(my_range,CDF, 'o--');
title('CDF - discrete uniform')
PMF = diff(CDF);
subplot(2,2,2);
stem(my_range(2:end), PMF, 'r');
title('PMF - discrete uniform')

%Discrete triangular
CF2 = CF1 * CF1;
my_range = 2*a-3:2*b+3;
CDF = Davies_inversion_discrete( CF2, my_range );
subplot(2,2,3);
stairs(my_range,CDF, 'o--');
title('CDF - discrete triangular')
PMF = diff(CDF);
subplot(2,2,4);
stem(my_range(2:end), PMF, 'r');
title('PMF - discrete triangular')