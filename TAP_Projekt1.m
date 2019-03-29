clear all;
close all;

%wpisanie danych
C=0.4;
a=8;
TCp=21;
THp=64;
TDp=30;
FCp=32;
FHp=20;
FDp=9;
tau=40;
tau_C=100;
hp=58.14;
Tp=36.43;

%% punkt 1
%% podpunkta)
%Symulacja obiektu - równania ró¿niczkowe\

%punkt pracy, skok tylko wartosci FH, skok tylko wartosci FC, skok obu wartosci 
FH=[0, 20;
    300, 22;
    900, 24];

FC=[0, 32;
    600, 34;
    900, 36];

FD = [0, 0];

x_0 = [54 35];

sim_time = 1200;
[t,x,y] = objectSimulation(FH,FC,FD,0,tau_C, 0, tau, THp, TCp, TDp, sim_time, x_0);

figure;
plot(t,x);
legend('h','T');
title('Przebieg zmiennych stanu (symulacja)');

figure;
plot(t,y);
legend('T_out')
title('Przebieg wyjscia obiektu (symulacja)');

%Obliczanie transmitancji
%Wyznaczanie wartoœci liczbowych macierzy A
a11=-(FHp+FCp+FDp)/(C*hp^2);
a12=((THp-Tp)*FHp+(TCp-Tp)*FCp+(TDp-Tp)*FDp)*(-2)/(C*hp^3);
a21=0;
a22=((FHp+FCp+FDp)/(hp^2)+a/(2*(sqrt(hp^3))))/(2*C);

A=[a11 a12; a21 a22];

%Wyznaczanie wartoœci liczbowych macierzy B
b11=1/(2*C*hp);
b12=1/(2*C*hp);
b21=THp-Tp;
b22=TCp-Tp;

B=[b11 b12; b21 b22];

%Wyznaczanie wartoœci liczbowych macierzy C
C=[0 1; 1 0];

%Wyznaczanie wartoœci macierzy D
D=[0 0; 0 0];

%Obliczenie transmitancji
syms s;
[G1]=C*(inv(s*eye(2)-A))*B+D;
% 
% %Uproszczenie równañ transmitancji
% G1=collect(G1);
% 
% 
% [numerator,denumerator]=numden(G1(1,1));
% G1n11=coeffs(numerator,s);
% G1dn11=coeffs(denumerator,s);
% 
% 
% [N, D]=numden(G1);
% %[N, D]=tfdata(G1(1,1));
