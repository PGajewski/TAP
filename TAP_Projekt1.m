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
Ts = 0.1;

%% punkt 1
%% podpunkt a)
%Symulacja obiektu - równania ró¿niczkowe\

%punkt pracy, skok tylko wartosci FH, skok tylko wartosci FC, skok obu wartosci 
FH=[0, 20;
    300, 30];

FC=[0,32];

FD = [0, 9];

TD = [0, 30];

sim_time = 1200;
[t1,x1,y1] = objectSimulation(FH,FC,FD,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);

figure(1);
plot(t1,x1);
legend('T','h');
title('Przebieg zmiennych stanu (symulacja)');

figure(2);
plot(t1,y1);
legend('T_{out}','h')
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

for k=1:4
    [num,den]=tfdata(G1(1,k));  %dodanie opoznienia do transmitancji GTi
    G1(1,k) = tf(num,den,'InputDelay',tau);
end

[num,den]=tfdata(G1(1,2));  %dodanie opoznienia do transmitancji GTFC
G1(1,2) = tf(num,den,'InputDelay',(tau+tau_C));
[num,den]=tfdata(G1(2,2));  %dodanie opoznienia do transmitancji GhFC
G1(2,2) = tf(num,den,'InputDelay',tau_C);

sys1 = ss(A,B,C,D,'InputDelay',[0,tau_C,0,0],'OutputDelay',[tau,0]);

FH=[0, 0;
    300, 10];

FC=[0, 0];

FD = [0, 0];

TD = [0, 0];

%przeprowadzenie symulacji
figure;
% %Dodanie opoznienia do transmitancji GxFC
% [num,den]=tfdata(G1);  %dodanie opoznienia do transmitancji GTi
% G = tf(num, den,'InputDelay',[0 tau_C, 0,0],'OutputDelay',[tau,0]);
% 
% %Przeprowadzenie symulacji
% sys_con = ss(A,B,C,D,'InputDelay',[0 tau_C, 0,0],'OutputDelay',[tau,0]);
% plot(t2,y2+ones(length(y2),2).*[Tp hp],t1',y1');



% 
%% punkt 1
%% podpunkt c)
% Transmitancja dyskretna
%G2 = c2d(G1, Ts, 'zoh');
%sys_discrete = c2d(sys_con,Ts,'zoh');
%[t3,x3,y3] = linearSimulation(FH,FC,FD,TD, sys_discrete, Ts, sim_time, [Tp, hp]);

%% Symulacje;
sim_time = 800;
FH=[0, 20;
    100 30];

FC=[0,32];

FD = [0, 9];

TD = [0, 30];

delta_FH=[0, 0;
          100, 10];

delta_FC=[0,0];

delta_FD = [0, 0];

delta_TD = [0, 0];

[t11,x11,y11] = objectSimulation(FH,FC,FD,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t12,x12,y12] = linearSimulation(delta_FH,delta_FC,delta_FD,delta_TD,sys1, Ts, sim_time, [0, 0]);

[t21,x21,y21] = objectSimulation(FH+[0 0; 0 10],FC,FD,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t22,x22,y22] = linearSimulation(delta_FH+[0 0; 0 10],delta_FC,delta_FD,delta_TD,sys1, Ts, sim_time, [0, 0]);

[t31,x31,y31] = objectSimulation(FH+[0 0; 0 20],FC,FD,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t32,x32,y32] = linearSimulation(delta_FH+[0 0; 0 20],delta_FC,delta_FD,delta_TD,sys1, Ts, sim_time, [0, 0]);


%% Rysowanie wykresów
% Porównanie zmiennych stanu.
figure(3);
plot(t11,y11(1,:),'r-',t12,y12(:,1)+ones(1,length(y12))*Tp,'b-',t21,y21(1,:),'r-',t22,y22(:,1)+ones(1,length(y22))*Tp,'b-',t31,y31(1,:),'r-',t32,y32(:,1)+ones(1,length(y32))*Tp,'b-');
legend('T nieliniowe','T zlinearyzowane');
title('Porownanie wyjœcia obiektu - temperatura');

