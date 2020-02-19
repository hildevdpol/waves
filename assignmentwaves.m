%% hilde van der pol 
%dec 2018

clear all
close all 
clc

%there are two strings, with variables: 
f = 12 ; 
w= 2*pi*f;
A = 3*10^-3; 
T = 1/f; 
v1 = 20 ; 
v2a = 0.3*v1 ; 
 
 
labda1 = v1*T; 
labda2 = v2a*T; 
 
k1 = 2*pi / labda1; 
k2 = 2*pi / labda2; 

r = (v2a-v1)/(v2a+v1); 
tr = (2*v2a)/(v2a+v1); 

t=0; 
 
  
%% Opdracht 1a
% Hilde van der Pol - 4663209
% plot one deviation for t= 0 and x=30/k to x = -30/k
 

x1 = (-30/k1:0.001:0); 
x2 = (0.001:0.001:30/k2) ;
 
 
for i = (1:1: length(x1))
        yin(i) = A*sin(k1*x1(i)-w*t);
        yr(i) = A*r*sin(-k1*x1(i)-w*t);
end
 
for i = ( 1:1: length(x2))
      
        ytr(i) = A*tr*sin(k2*x2(i)-w*t);
end
 
 ytotb = yin+yr;
   

figure(1)
hold on 
xlabel('afstand x [m]')
ylabel('uitwijking y [m]')
title ('koord op t=0')

plot (x1,ytotb,'k')
plot (x2,ytr, 'k')
 
 
%% Opdracht 2
% Hilde van der Pol - 4663209
%do the same for t= 1/4f 

t1b= 1/(4*f); 
 
 
for i = (1:1: length(x1))
        yinb(i) = A*sin(k1*x1(i)-w*t1b);
        yrb(i) = A*r*sin(-k1*x1(i)-w*t1b);
end
 
for i = ( 1:1: length(x2))
      
        ytrb(i) = A*tr*sin(k2*x2(i)-w*t1b);
end
 
ytotb = yinb+yrb; 

figure(2)
hold on 
xlabel('afstand x [m]')
ylabel('uitwijking y [m]')
title('Koord op t= 1/4f')

plot (x1,ytotb, 'k')
plot (x2,ytrb, 'k')
 
%% Opdracht 3
% Hilde van der Pol, 4663209
%plot for x=0 the total deviation of the string from t=-5/f to t-5/f
 
x3= 0; 
t3= (-5/f: 0.001: 5/f); 
 for i= (1: 1: length(t3))
    yin3(i)= A*sin(k1*x3-w*t3(i));
yr3(i)= A*r*sin(-k1*x3-w*t3(i));
end

ytot3 = yin3 + yr3;

figure(3)
hold on 
xlabel('tijd t [s]')
ylabel('uitwijking y [m]')
title('Koord op x=0')

plot (t3, ytot3, 'k')
 
 
 
%% Opdracht 4
% Hilde van der Pol
% do the same for:x= n/k1 with n = -4, -3, -2, -1 and for x = n/k2 with n = 1,2,3,4   

n1 = [-4, -3, -2, -1];
x41 = n1/k1;

figure (4)
xlabel('Uitwijking x [m]')
ylabel('Maximale amplitude [m]')
title('Totale uitwijking koord op verschillende posities x')   

hold on 
for ii= (1:1: length(x41))
    for jj = (1:1: length(t3))
        yin4(ii,jj) = A*sin(k1*x41(ii)-w*t3(jj));
        yr4(ii,jj) = A*r*sin(-k1*x41(ii)-w*t3(jj));
    end
    hold on
   
    ytot4= yin4 + yr4;
    plot(t3,ytot4);
end
 
n2 = [1, 2, 3, 4];
x42 = n2/k2; 


for ii = (1: 1: length(x42))
    for jj = (1: 1: length(t3))
        yt4(ii,jj) = A*tr*sin((k2*x42(ii))-w*t3(jj));
    end
    hold on
    plot (t3,yt4);
end
 

 
%% Opgave 5
% Hilde van der Pol, 4663209
%plot the amplitude of the total wave from x= -10/k1 to x= 10/k2 for: 
% v2= 0.3v1, v2=v1, v2=3v1
 
t5 = (-1/f: 0.001: 1/f);
  
v2b = v1 ; 
v2c = 3*v1 ;  

labda2b = v2b*T; 
labda2c = v2c*T; 
k2b = 2*pi / labda2b; 
k2c = 2*pi / labda2c; 
rb = (v2b-v1)/(v2b+v1); 
rc = (v2c-v1)/(v2c+v1); 
trb = (2*v2b)/(v2b+v1); 
trc = (2*v2c)/(v2c+v1); 

x5a = ((-10/k1):0.01:0);
x5b = (0.00001:0.01:(10/k2));
x5bb = (0.00001:0.01:(10/k2b));
x5bc = (0.00001:0.01:(10/k2c));



figure(5)
hold on
%Voor v2 = 0.3v1
% Hilde van der Pol, 4663209
for jj = (1: 1: length(t5))
    for ii = (1: 1: length(x5a))
        yin5(jj,ii) = A*sin(k1*x5a(ii)-w*t5(jj));
        yr5(jj,ii) = A*r*sin(-k1*x5a(ii)-w*t5(jj));
        ytot5(jj,ii) = yin5(jj,ii)+yr5(jj,ii);
        maxytot5 = max(ytot5);
    end
    
    hold on
end
 plot (x5a, maxytot5, 'b')
 
for jj = (1: 1: length(t5))
    for ii = (1: 1: length(x5b))
       yt5(jj,ii)= A*tr*sin((k2*x5b(ii))-w*t5(jj));
        maxyt5 = max(yt5);
    end
   
end
 plot (x5b, maxyt5, 'b') 

%Voor v2 = v1
% Hilde van der Pol, 4663209
for jj = (1: 1: length(t5))
    for ii = (1: 1: length(x5a))
        yin6(jj,ii) = A*sin(k1*x5a(ii)-w*t5(jj));
        yr6(jj,ii) = A*rb*sin(-k1*x5a(ii)-w*t5(jj));
        ytot6(jj,ii) = yin6(jj,ii)+yr6(jj,ii);
        maxytot6 = max(ytot6);
    end
   
end
 plot (x5a, maxytot6, 'g')
 
hold on
for jj = (1: 1: length(t5))
    for ii = (1: 1: length(x5bb))
       yt6(jj,ii)= A*trb*sin((k2b*x5bb(ii))-w*t5(jj));
        maxyt6 = max(yt6);
    end
   
end
 plot (x5bb, maxyt6, 'g') 

% Voor v2 = 3v1
% Hilde van der Pol, 4663209
for jj = (1: 1: length(t5))
    for ii = (1: 1: length(x5a))
        yin7(jj,ii) = A*sin(k1*x5a(ii)-w*t5(jj));
        yr7(jj,ii) = A*rc*sin(-k1*x5a(ii)-w*t5(jj));
        ytot7(jj,ii) = yin7(jj,ii)+yr7(jj,ii);
        maxytot7 = max(ytot7);
    end
   
end
 plot (x5a, maxytot7, 'r')
 
hold on
for jj = (1: 1: length(t5))
    for ii = (1: 1: length(x5bc))
       yt7(jj,ii)= A*trc*sin((k2c*x5bc(ii))-w*t5(jj));
        maxyt7 = max(yt7);
    end
    
end
plot (x5bc, maxyt7, 'r')

title 'Amplitude totale golf voor verschillende v'
xlabel 'Verplaatsing x [m]'
ylabel 'maximale amplitude [m]'