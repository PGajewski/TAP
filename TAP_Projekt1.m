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
Ts = 10;

%% punkt 1
%% podpunkt a)
%Symulacja obiektu - równania ró¿niczkowe\

%punkt pracy, skok tylko wartosci FH, skok tylko wartosci FC, skok obu wartosci 
FH=[0, 20;
    300, 22;
    900, 24];

FC=[0, 32;
    600, 34;
    900, 36];

FD = [0, 9];

TD = [0, 30];

x_0 = [54 35];

sim_time = 1200;
[t1,x1,y1] = objectSimulation(FH,FC,FD,TD,0,tau_C, 0, tau, THp, TCp, TDp, sim_time, x_0);

figure;
plot(t1,x1);
legend('h','T');
title('Przebieg zmiennych stanu (symulacja)');

figure;
plot(t1,y1);
legend('T_out')
title('Przebieg wyjscia obiektu (symulacja)');

%Obliczanie transmitancji
%Wyznaczanie wartoœci liczbowych macierzy A
a11=-(FHp+FCp+FDp)/(C*hp^2);
a12=(-2)*((THp-Tp)*FHp+(TCp-Tp)*FCp+(TDp-Tp)*FDp)/(C*hp^3);
a21=0;
a22=((FHp+FCp+FDp)/(hp^2)+a/(2*(sqrt(hp^3))))/(-2*C);

A=[a11 a12; a21 a22];

%Wyznaczanie wartoœci liczbowych macierzy B
b11=1/(C*hp^2)*(THp-Tp);
b12=1/(C*hp^2)*(TCp-Tp);
b13=1/(C*hp^2)*(TDp-Tp);
b14=1/(C*hp^2)*FDp;
b21=1/(2*C*hp);
b22=1/(2*C*hp);
b23=1/(2*C*hp);
b24=0;

B=[b11 b12 b13 b14; b21 b22 b23 b24];

%Wyznaczanie wartoœci liczbowych macierzy C
C=[1 0; 0 1];

%Wyznaczanie wartoœci macierzy D
D=[0 0 0 0; 0 0 0 0];

%Obliczenie transmitancji
syms s;
tf('s');
G1=C*(inv(s*eye(2)-A))*B+D;
% 
% %Uproszczenie równañ transmitancji
G1=collect(G1);
G11=G1(1,1);
G12=G1(1,2);
G13=G1(1,3);
G14=G1(1,4);

G21=G1(2,1);
G22=G1(2,2);
G23=G1(2,3);
G24=G1(2,4);
clear G1;

 %ustalenie s jako zmiennej transformaty
s = tf('s');
 
%Przekszta³cenie wielomianów na w³aœciwe transformaty
G=char(G11);
eval(['G1(1,1)=',G]);
G=char(G12);
eval(['G1(1,2)=',G]);
G=char(G13);
eval(['G1(1,3)=',G]);
G=char(G14);
eval(['G1(1,4)=',G]);

G=char(G21);
eval(['G1(2,1)=',G]);
G=char(G22);
eval(['G1(2,2)=',G]);
G=char(G23);
eval(['G1(2,3)=',G]);
G=char(G24);
eval(['G1(2,4)=',G]);

%Dodanie opóŸnienia na wyjœciu.
for k=1:4
    [num,den]=tfdata(G1(1,k));  %dodanie opoznienia do transmitancji GTi
    G1(1,k) = tf(num,den,'OutputDelay',tau);
end

%Dodanie opoznienia do transmitancji GxFC
for k=1:2
    G1(k,2) = tf(num,den,'InputDelay',(tau_C));
end

%Przeprowadzenie symulacji
[t2,x2,y2] = linearSimulation(FH, FC, FD, TD, G1, Ts, sim_time, x_0);

%% punkt 1
%% podpunkt c)
% Transmitancja dyskretna
G2 = c2d(G1, Ts, 'zoh');

[t3,x3,y3] = linearSimulation(FH,FC,FD,TD, G2, Ts, sim_time, x_0)

%% Rysowanie wykresów
% Porównanie zmiennych stanu.
figure(1);
plot(t1,x1(1,:),t2,x2(:,1)');
figure(2);
plot(t1,x1(2,:),t2,x2(:,2)');
figure =(3);
plot(t1,y1,t2',y2');

