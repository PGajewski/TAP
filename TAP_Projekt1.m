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
%Symulacja obiektu - rownania rowniczkowe.

%punkt pracy, skok tylko wartosci FH, skok tylko wartosci FC, skok obu wartosci 
FH=[0, 20;
    300, 30];

FC=[0,32;
    400, 35];

FD = [0, 9;
      700, 13];

TD = [0, 30;
      900 40];

sim_time = 1200;
[t1,x1,y1] = objectSimulation(FH,FC,FD,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);

figure(1);
plot(t1,x1);
legend('T','h');
xlabel('Czas [s]');
ylabel('Stany wewnetrzne');
title('Przebieg zmiennych stanu (symulacja)');

figure(2);
plot(t1,y1);
legend('T_{out}','h')
xlabel('czas [s]');
ylabel('Wyjscia obiektu');
title('Przebieg wyjscia obiektu (symulacja)');

%Obliczanie transmitancji
%Wyznaczanie wartosci liczbowych macierzy A
a11=-(FHp+FCp+FDp)/(C*hp^2);
a12=(-2)*((THp-Tp)*FHp+(TCp-Tp)*FCp+(TDp-Tp)*FDp)/(C*hp^3);
a21=0;
a22=((FHp+FCp+FDp)/(hp^2)-a/(2*(sqrt(hp^3))))/(-2*C);

A=[a11 a12; a21 a22];

%Wyznaczanie wartosci liczbowych macierzy B
b11=1/(C*hp^2)*(THp-Tp);
b12=1/(C*hp^2)*(TCp-Tp);
b13=1/(C*hp^2)*(TDp-Tp);
b14=1/(C*hp^2)*FDp;
b21=1/(2*C*hp);
b22=1/(2*C*hp);
b23=1/(2*C*hp);
b24=0;

B=[b11 b12 b13 b14; b21 b22 b23 b24];

%Wyznaczanie wartosci liczbowych macierzy C
C=[1 0; 0 1];

%Wyznaczanie wartosci macierzy D
D=[0 0 0 0; 0 0 0 0];

%Obliczenie transmitancji
syms s;
tf('s');
G1=C*(inv(s*eye(2)-A))*B+D;
% 
% %Uproszczenie rownan transmitancji
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
 
%Przeksztalcenie wielomianow na wlasciwe transformaty
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

[num,den]=tfdata(G1);  %dodanie opoznienia do transmitancji GTFC
%G1(1,2) = tf(num,den,'InputDelay',(tau+tau_C));
%[num,den]=tfdata(G1(2,2));  %dodanie opoznienia do transmitancji GhFC
%G1(2,2) = tf(num,den,'InputDelay',tau_C);

G1 = tf(num,den,'InputDelay',[0,tau_C,0,0],'OutputDelay',[tau,0]);

sys_con = ss(A,B,C,D,'InputDelay',[0,tau_C,0,0],'OutputDelay',[tau,0]);

%Uklady dyskretne.
Ts1 = 10;
Ts2 = 20;
Ts3 = 20;

sys1 = c2d(sys_con,Ts1);
sys2 = c2d(sys_con,Ts2);
sys3 = c2d(sys_con,Ts3);

%% punkt 1
%% podpunkt c)
% Transmitancja dyskretna
tf1 = c2d(G1,Ts1);
tf2 = c2d(G1,Ts2);
tf3 = c2d(G1,Ts3);


%% Symulacje;
sim_time = 800;
sim_discrete_time = 800;
sim_step = [0 0; 0 5];
sim_start = 30;
FH=[0, FHp;
    100, 20];

FC = [0,FCp;
      100,32];

FD = [0,FDp;
      100,9];

TD = [0,TDp;
      100,30];

delta_FH=[0,0;
          100,0];

delta_FC = [0,0;
            100,0];

delta_FD = [0,0;
            100,0];

delta_TD = [0,0;
            100,0];
%Skok FH
%Symulacja przy skokach dodatnich
[t11,x11,y11] = objectSimulation(FH+sim_step,FC,FD,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t12,x12,y12] = linearSimulation(delta_FH+sim_step,delta_FC,delta_FD,delta_TD,sys_con, Ts, sim_time, [0, 0]);

[t21,x21,y21] = objectSimulation(FH+2*sim_step,FC,FD,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t22,x22,y22] = linearSimulation(delta_FH+2*sim_step,delta_FC,delta_FD,delta_TD,sys_con, Ts, sim_time, [0, 0]);

[t31,x31,y31] = objectSimulation(FH+3*sim_step,FC,FD,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t32,x32,y32] = linearSimulation(delta_FH+3*sim_step,delta_FC,delta_FD,delta_TD,sys_con, Ts, sim_time, [0, 0]);

%Symulacja przy skokach ujemnych
[t41,x41,y41] = objectSimulation(FH-sim_step,FC,FD,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t42,x42,y42] = linearSimulation(delta_FH-sim_step,delta_FC,delta_FD,delta_TD,sys_con, Ts, sim_time, [0, 0]);

[t51,x51,y51] = objectSimulation(FH-2*sim_step,FC,FD,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t52,x52,y52] = linearSimulation(delta_FH-2*sim_step,delta_FC,delta_FD,delta_TD,sys_con, Ts, sim_time, [0, 0]);

[t61,x61,y61] = objectSimulation(FH-3*sim_step,FC,FD,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t62,x62,y62] = linearSimulation(delta_FH-3*sim_step,delta_FC,delta_FD,delta_TD,sys_con, Ts, sim_time, [0, 0]);

%Skok FC
%Symulacja przy skokach dodatnich
[t71,x71,y71] = objectSimulation(FH,FC+sim_step,FD,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t72,x72,y72] = linearSimulation(delta_FH,delta_FC+sim_step,delta_FD,delta_TD,sys_con, Ts, sim_time, [0, 0]);

[t81,x81,y81] = objectSimulation(FH,FC+2*sim_step,FD,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t82,x82,y82] = linearSimulation(delta_FH,delta_FC+2*sim_step,delta_FD,delta_TD,sys_con, Ts, sim_time, [0, 0]);

[t91,x91,y91] = objectSimulation(FH,FC+3*sim_step,FD,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t92,x92,y92] = linearSimulation(delta_FH,delta_FC+3*sim_step,delta_FD,delta_TD,sys_con, Ts, sim_time, [0, 0]);

%Symulacja przy skokach ujemnych
[t101,x101,y101] = objectSimulation(FH,FC-sim_step,FD,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t102,x102,y102] = linearSimulation(delta_FH,delta_FC-sim_step,delta_FD,delta_TD,sys_con, Ts, sim_time, [0, 0]);

[t111,x111,y111] = objectSimulation(FH,FC-2*sim_step,FD,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t112,x112,y112] = linearSimulation(delta_FH,delta_FC-2*sim_step,delta_FD,delta_TD,sys_con, Ts, sim_time, [0, 0]);

[t121,x121,y121] = objectSimulation(FH,FC-3*sim_step,FD,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t122,x122,y122] = linearSimulation(delta_FH,delta_FC-3*sim_step,delta_FD,delta_TD,sys_con, Ts, sim_time, [0, 0]);

%Skok FD
%Symulacja przy skokach dodatnich
[t131,x131,y131] = objectSimulation(FH,FC,FD+0.2*sim_step,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t132,x132,y132] = linearSimulation(delta_FH,delta_FC,delta_FD+0.2*sim_step,delta_TD,sys_con, Ts, sim_time, [0, 0]);

[t141,x141,y141] = objectSimulation(FH,FC,FD+0.4*sim_step,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t142,x142,y142] = linearSimulation(delta_FH,delta_FC,delta_FD+0.4*sim_step,delta_TD,sys_con, Ts, sim_time, [0, 0]);

[t151,x151,y151] = objectSimulation(FH,FC,FD+0.6*sim_step,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t152,x152,y152] = linearSimulation(delta_FH,delta_FC,delta_FD+0.6*sim_step,delta_TD,sys_con, Ts, sim_time, [0, 0]);

%Symulacja przy skokach ujemnych
[t161,x161,y161] = objectSimulation(FH,FC,FD-0.2*sim_step,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t162,x162,y162] = linearSimulation(delta_FH,delta_FC,delta_FD-0.2*sim_step,delta_TD,sys_con, Ts, sim_time, [0, 0]);

[t171,x171,y171] = objectSimulation(FH,FC,FD-0.4*sim_step,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t172,x172,y172] = linearSimulation(delta_FH,delta_FC,delta_FD-0.4*sim_step,delta_TD,sys_con, Ts, sim_time, [0, 0]);

[t181,x181,y181] = objectSimulation(FH,FC,FD-0.6*sim_step,TD,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t182,x182,y182] = linearSimulation(delta_FH,delta_FC,delta_FD-0.6*sim_step,delta_TD,sys_con, Ts, sim_time, [0, 0]);

%Skok TD
%Symulacja przy skokach dodatnich
[t191,x191,y191] = objectSimulation(FH,FC,FD,TD+sim_step,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t192,x192,y192] = linearSimulation(delta_FH,delta_FC,delta_FD,delta_TD+sim_step,sys_con, Ts, sim_time, [0, 0]);

[t201,x201,y201] = objectSimulation(FH,FC,FD,TD+2*sim_step,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t202,x202,y202] = linearSimulation(delta_FH,delta_FC,delta_FD,delta_TD+2*sim_step,sys_con, Ts, sim_time, [0, 0]);

[t211,x211,y211] = objectSimulation(FH,FC,FD,TD+3*sim_step,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t212,x212,y212] = linearSimulation(delta_FH,delta_FC,delta_FD,delta_TD+3*sim_step,sys_con, Ts, sim_time, [0, 0]);

%Symulacja przy skokach ujemnych
[t221,x221,y221] = objectSimulation(FH,FC,FD,TD-sim_step,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t222,x222,y222] = linearSimulation(delta_FH,delta_FC,delta_FD,delta_TD-sim_step,sys_con, Ts, sim_time, [0, 0]);

[t231,x231,y231] = objectSimulation(FH,FC,FD,TD-2*sim_step,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t232,x232,y232] = linearSimulation(delta_FH,delta_FC,delta_FD,delta_TD-2*sim_step,sys_con, Ts, sim_time, [0, 0]);

[t241,x241,y241] = objectSimulation(FH,FC,FD,TD-3*sim_step,0,tau_C, 0, tau, THp, TCp, sim_time, [Tp, hp]);
[t242,x242,y242] = linearSimulation(delta_FH,delta_FC,delta_FD,delta_TD-3*sim_step,sys_con, Ts, sim_time, [0, 0]);

%Symulacje dyskretne
%Ts = 10
[t251,x251,y251] = linearSimulation(delta_FH+sim_step,delta_FC,delta_FD,delta_TD,sys1, Ts1, sim_discrete_time, [0, 0]);
[t252,x252,y252] = linearSimulation(delta_FH,delta_FC+sim_step,delta_FD,delta_TD,sys1, Ts1, sim_discrete_time, [0, 0]);
[t253,x253,y253] = linearSimulation(delta_FH,delta_FC,delta_FD+0.2*sim_step,delta_TD,sys1, Ts1, sim_discrete_time, [0, 0]);
[t254,x254,y254] = linearSimulation(delta_FH,delta_FC,delta_FD,delta_TD+sim_step,sys1, Ts1, sim_discrete_time, [0, 0]);

%Ts = 20
[t261,x261,y261] = linearSimulation(delta_FH+sim_step,delta_FC,delta_FD,delta_TD,sys2, Ts2, sim_discrete_time, [0, 0]);
[t262,x262,y262] = linearSimulation(delta_FH,delta_FC+sim_step,delta_FD,delta_TD,sys2, Ts2, sim_discrete_time, [0, 0]);
[t263,x263,y263] = linearSimulation(delta_FH,delta_FC,delta_FD+0.2*sim_step,delta_TD,sys2, Ts2, sim_discrete_time, [0, 0]);
[t264,x264,y264] = linearSimulation(delta_FH,delta_FC,delta_FD,delta_TD+sim_step,sys2, Ts2, sim_discrete_time, [0, 0]);

%Ts = 25
[t271,x271,y271] = linearSimulation(delta_FH+sim_step,delta_FC,delta_FD,delta_TD,sys3, Ts3, sim_discrete_time, [0, 0]);
[t272,x272,y272] = linearSimulation(delta_FH,delta_FC+sim_step,delta_FD,delta_TD,sys3, Ts3, sim_discrete_time, [0, 0]);
[t273,x273,y273] = linearSimulation(delta_FH,delta_FC,delta_FD+0.2*sim_step,delta_TD,sys3, Ts3, sim_discrete_time, [0, 0]);
[t274,x274,y274] = linearSimulation(delta_FH,delta_FC,delta_FD,delta_TD+sim_step,sys3, Ts3, sim_discrete_time, [0, 0]);

%Symulacje transmitancji 
%Model ciagly
[t281,x281,y281] = linearSimulation(delta_FH+sim_step,delta_FC,delta_FD,delta_TD,G1, 0.1, sim_time, [0, 0]);
[t282,x282,y282] = linearSimulation(delta_FH,delta_FC+sim_step,delta_FD,delta_TD,G1, 0.1, sim_time, [0, 0]);
[t283,x283,y283] = linearSimulation(delta_FH,delta_FC,delta_FD+0.2*sim_step,delta_TD,G1, 0.1, sim_time, [0, 0]);
[t284,x284,y284] = linearSimulation(delta_FH,delta_FC,delta_FD,delta_TD+sim_step,G1, 0.1, sim_time, [0, 0]);

%Model dyskretny, Ts = 10
[t291,x291,y291] = linearSimulation(delta_FH+sim_step,delta_FC,delta_FD,delta_TD,tf1, Ts1, sim_discrete_time, [0, 0]);
[t292,x292,y292] = linearSimulation(delta_FH,delta_FC+sim_step,delta_FD,delta_TD,tf1, Ts1, sim_discrete_time, [0, 0]);
[t293,x293,y293] = linearSimulation(delta_FH,delta_FC,delta_FD+0.2*sim_step,delta_TD,tf1, Ts1, sim_discrete_time, [0, 0]);
[t294,x294,y294] = linearSimulation(delta_FH,delta_FC,delta_FD,delta_TD+sim_step,tf1, Ts1, sim_discrete_time, [0, 0]);


%Model dyskretny, Ts = 20
[t301,x301,y301] = linearSimulation(delta_FH+sim_step,delta_FC,delta_FD,delta_TD,tf2, Ts2, sim_discrete_time, [0, 0]);
[t302,x302,y302] = linearSimulation(delta_FH,delta_FC+sim_step,delta_FD,delta_TD,tf2, Ts2, sim_discrete_time, [0, 0]);
[t303,x303,y303] = linearSimulation(delta_FH,delta_FC,delta_FD+0.2*sim_step,delta_TD,tf2, Ts2, sim_discrete_time, [0, 0]);
[t304,x304,y304] = linearSimulation(delta_FH,delta_FC,delta_FD,delta_TD+sim_step,tf2, Ts2, sim_discrete_time, [0, 0]);

%Model dyskretny, Ts = 25
[t311,x311,y311] = linearSimulation(delta_FH+sim_step,delta_FC,delta_FD,delta_TD,tf3, Ts3, sim_discrete_time, [0, 0]);
[t312,x312,y312] = linearSimulation(delta_FH,delta_FC+sim_step,delta_FD,delta_TD,tf3, Ts3, sim_discrete_time, [0, 0]);
[t313,x313,y313] = linearSimulation(delta_FH,delta_FC,delta_FD+0.2*sim_step,delta_TD,tf3, Ts3, sim_discrete_time, [0, 0]);
[t314,x314,y314] = linearSimulation(delta_FH,delta_FC,delta_FD,delta_TD+sim_step,tf3, Ts3, sim_discrete_time, [0, 0]);
%%
% 
% $$e^{\pi i} + 1 = 0$$
% 


%% Rysowanie wykresow
% Porownanie wyjsc obiektu nieliniowego oraz modelu zlinearyzowanego, przy skokach dodatnich
% Skok FH
% Skok dodatni
figure(3);
plot(t11(1,:),y11(1,:),'r-',t12(:,1),y12(:,1)+(ones(1,length(y12))*Tp)','b-',t21(1,:),y21(1,:),'r-',t22(:,1),y22(:,1)+(ones(1,length(y22))*Tp)','b-',t31(1,:),y31(1,:),'r-',t32(:,1),y32(:,1)+(ones(1,length(y32))*Tp)','b-');
legend('obiekt nieliniowy','model liniowy');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Skok dodatni FH - temperatura');

figure(4);
plot(t11(1,:),y11(2,:),'r-',t12(:,1),y12(:,2)+(ones(1,length(y12))*hp)','b-',t21(1,:),y21(2,:),'r-',t22(:,1),y22(:,2)+(ones(1,length(y22))*hp)','b-',t31(1,:),y31(2,:),'r-',t32(:,1),y32(:,2)+(ones(1,length(y32))*hp)','b-');
legend('obiekt nieliniowy','model liniowy');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Skok dodatni FH - wysokosc');

%Skok ujemny
figure(5);
plot(t41(1,:),y41(1,:),'r-',t42(:,1),y42(:,1)+(ones(1,length(y42))*Tp)','b-',t51(1,:),y51(1,:),'r-',t52(:,1),y52(:,1)+(ones(1,length(y52))*Tp)','b-',t61(1,:),y61(1,:),'r-',t62(:,1),y62(:,1)+(ones(1,length(y62))*Tp)','b-');
legend('obiekt nieliniowy','model liniowy');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Skok ujemny FH - temperatura');

figure(6);
plot(t41(1,:),y41(2,:),'r-',t42(:,1),y42(:,2)+(ones(1,length(y42))*hp)','b-',t51(1,:),y51(2,:),'r-',t52(:,1),y52(:,2)+(ones(1,length(y52))*hp)','b-',t61(1,:),y61(2,:),'r-',t62(:,1),y62(:,2)+(ones(1,length(y62))*hp)','b-');
legend('obiekt nieliniowy','model liniowy');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Skok ujemny FH  - wysokosc');

% Skok FC
% Skok dodatni
figure(7);
plot(t71(1,:),y71(1,:),'r-',t72(:,1),y72(:,1)+(ones(1,length(y72))*Tp)','b-',t81(1,:),y81(1,:),'r-',t82(:,1),y82(:,1)+(ones(1,length(y82))*Tp)','b-',t91(1,:),y91(1,:),'r-',t92(:,1),y92(:,1)+(ones(1,length(y92))*Tp)','b-');
legend('obiekt nieliniowy','model liniowy');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Skok dodatni FC  - temperatura');

figure(8);
plot(t71(1,:),y71(2,:),'r-',t72(:,1),y72(:,2)+(ones(1,length(y72))*hp)','b-',t81(1,:),y81(2,:),'r-',t82(:,1),y82(:,2)+(ones(1,length(y82))*hp)','b-',t91(1,:),y91(2,:),'r-',t92(:,1),y92(:,2)+(ones(1,length(y92))*hp)','b-');
legend('obiekt nieliniowy','model liniowy');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Skok dodatni FC - wysokosc');

%Skok ujemny
figure(9);
plot(t101(1,:),y101(1,:),'r-',t102(:,1),y102(:,1)+(ones(1,length(y102))*Tp)','b-',t111(1,:),y111(1,:),'r-',t112(:,1),y112(:,1)+(ones(1,length(y112))*Tp)','b-',t121(1,:),y121(1,:),'r-',t122(:,1),y122(:,1)+(ones(1,length(y122))*Tp)','b-');
legend('obiekt nieliniowy','model liniowy');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Skok ujemny FC - temperatura');

figure(10);
plot(t101(1,:),y101(2,:),'r-',t102(:,1),y102(:,2)+(ones(1,length(y102))*hp)','b-',t111(1,:),y111(2,:),'r-',t112(:,1),y112(:,2)+(ones(1,length(y112))*hp)','b-',t121(1,:),y121(2,:),'r-',t122(:,1),y122(:,2)+(ones(1,length(y122))*hp)','b-');
legend('obiekt nieliniowy','model liniowy');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Skok ujemny FC - wysokosc');

% Skok FD
% Skok dodatni
figure(11);
plot(t131(1,:),y131(1,:),'r-',t132(:,1),y132(:,1)+(ones(1,length(y132))*Tp)','b-',t141(1,:),y141(1,:),'r-',t142(:,1),y142(:,1)+(ones(1,length(y142))*Tp)','b-',t151(1,:),y151(1,:),'r-',t152(:,1),y152(:,1)+(ones(1,length(y152))*Tp)','b-');
legend('obiekt nieliniowy','model liniowy');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Skok dodatni FD  - temperatura');

figure(12);
plot(t131(1,:),y131(2,:),'r-',t132(:,1),y132(:,2)+(ones(1,length(y132))*hp)','b-',t141(1,:),y141(2,:),'r-',t142(:,1),y142(:,2)+(ones(1,length(y142))*hp)','b-',t151(1,:),y151(2,:),'r-',t152(:,1),y152(:,2)+(ones(1,length(y152))*hp)','b-');
legend('obiekt nieliniowy','model liniowy');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Skok dodatni FD - wysokosc');

%Skok ujemny
figure(13);
plot(t161(1,:),y161(1,:),'r-',t162(:,1),y162(:,1)+(ones(1,length(y162))*Tp)','b-',t171(1,:),y171(1,:),'r-',t172(:,1),y172(:,1)+(ones(1,length(y172))*Tp)','b-',t181(1,:),y181(1,:),'r-',t182(:,1),y182(:,1)+(ones(1,length(y182))*Tp)','b-');
legend('obiekt nieliniowy','model liniowy');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Skok ujemny FD - temperatura');

figure(14);
plot(t161(1,:),y161(2,:),'r-',t162(:,1),y162(:,2)+(ones(1,length(y162))*hp)','b-',t171(1,:),y171(2,:),'r-',t172(:,1),y172(:,2)+(ones(1,length(y172))*hp)','b-',t181(1,:),y181(2,:),'r-',t182(:,1),y182(:,2)+(ones(1,length(y182))*hp)','b-');
legend('obiekt nieliniowy','model liniowy');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Skok ujemny FD - wysokosc');

% Skok TD
% Skok dodatni
figure(15);
plot(t191(1,:),y191(1,:),'r-',t192(:,1),y192(:,1)+(ones(1,length(y192))*Tp)','b-',t201(1,:),y201(1,:),'r-',t202(:,1),y202(:,1)+(ones(1,length(y202))*Tp)','b-',t211(1,:),y211(1,:),'r-',t212(:,1),y212(:,1)+(ones(1,length(y212))*Tp)','b-');
legend('obiekt nieliniowy','model liniowy');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Skok dodatni TD  - temperatura');

figure(16);
plot(t191(1,:),y191(2,:),'r-',t192(:,1),y192(:,2)+(ones(1,length(y192))*hp)','b-',t201(1,:),y201(2,:),'r-',t202(:,1),y202(:,2)+(ones(1,length(y202))*hp)','b-',t211(1,:),y211(2,:),'r-',t212(:,1),y212(:,2)+(ones(1,length(y212))*hp)','b-');
legend('obiekt nieliniowy','model liniowy');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Skok dodatni TD - wysokosc');

%Skok ujemny
figure(17);
plot(t221(1,:),y221(1,:),'r-',t222(:,1),y222(:,1)+(ones(1,length(y222))*Tp)','b-',t231(1,:),y231(1,:),'r-',t232(:,1),y232(:,1)+(ones(1,length(y232))*Tp)','b-',t241(1,:),y241(1,:),'r-',t242(:,1),y242(:,1)+(ones(1,length(y242))*Tp)','b-');
legend('obiekt nieliniowy','model liniowy');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Skok ujemny TD - temperatura');

figure(18);
plot(t221(1,:),y221(2,:),'r-',t222(:,1),y222(:,2)+(ones(1,length(y222))*hp)','b-',t231(1,:),y231(2,:),'r-',t232(:,1),y232(:,2)+(ones(1,length(y232))*hp)','b-',t241(1,:),y241(2,:),'r-',t242(:,1),y242(:,2)+(ones(1,length(y242))*hp)','b-');
legend('obiekt nieliniowy','model liniowy');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Skok ujemny TD - wysokosc');

%Porownanie bledow dyskretyzacji w modelach dyskretnych o roznych okresach probkowania
%Skok FH
figure(19)
stairs(t251(:,1),y251(:,1)+(ones(1,length(y251))*Tp)');
hold on;
stairs(t261(:,1),y261(:,1)+(ones(1,length(y261))*Tp)');
hold on;
stairs(t271(:,1),y271(:,1)+(ones(1,length(y271))*Tp)');
hold on;
plot(t12(:,1),y12(:,1)+(ones(1,length(y12))*Tp)','m-');
legend('Ts = 10','Ts = 20','Ts = 25','Uklad ciagly');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Skok dodatni FH dla roznych ukladow dyskretnych - temperatura');

figure(20)
stairs(t251(:,1),y251(:,2)+(ones(1,length(y251))*hp)');
hold on;
stairs(t261(:,1),y261(:,2)+(ones(1,length(y261))*hp)');
hold on;
stairs(t271(:,1),y271(:,2)+(ones(1,length(y271))*hp)');
hold on;
plot(t12(:,1),y12(:,2)+(ones(1,length(y12))*hp)','m-');
legend('Ts = 10','Ts = 20','Ts = 25','Uklad ciagly');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Skok dodatni FH dla roznych ukladow dyskretnych - wysokosc');

%Skok FC
figure(21)
stairs(t252(:,1),y252(:,1)+(ones(1,length(y252))*Tp)');
hold on;
stairs(t262(:,1),y262(:,1)+(ones(1,length(y262))*Tp)');
hold on;
stairs(t272(:,1),y272(:,1)+(ones(1,length(y272))*Tp)');
hold on;
plot(t72(:,1),y72(:,1)+(ones(1,length(y72))*Tp)','m-');
legend('Ts = 10','Ts = 20','Ts = 25','Uklad ciagly');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Skok dodatni FC dla roznych ukladow dyskretnych - temperatura');

figure(22)
stairs(t252(:,1),y251(:,2)+(ones(1,length(y252))*hp)');
hold on;
stairs(t262(:,1),y261(:,2)+(ones(1,length(y262))*hp)');
hold on;
stairs(t272(:,1),y271(:,2)+(ones(1,length(y272))*hp)');
hold on;
plot(t72(:,1),y12(:,2)+(ones(1,length(y72))*hp)','m-');
legend('Ts = 10','Ts = 20','Ts = 25','Uklad ciagly');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Skok dodatni FC dla roznych ukladow dyskretnych - wysokosc');

%Skok FD
figure(23)
stairs(t253(:,1),y253(:,1)+(ones(1,length(y253))*Tp)');
hold on;
stairs(t263(:,1),y263(:,1)+(ones(1,length(y263))*Tp)');
hold on;
stairs(t273(:,1),y273(:,1)+(ones(1,length(y273))*Tp)');
hold on;
plot(t132(:,1),y132(:,1)+(ones(1,length(y132))*Tp)','m-');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');legend('Ts = 10','Ts = 20','Ts = 25','Uklad ciagly');
title('Skok dodatni FD dla roznych ukladow dyskretnych - temperatura');

figure(24)
stairs(t253(:,1),y253(:,1)+(ones(1,length(y253))*hp)');
hold on;
stairs(t263(:,1),y263(:,1)+(ones(1,length(y263))*hp)');
hold on;
stairs(t273(:,1),y273(:,1)+(ones(1,length(y273))*hp)');
hold on;
plot(t132(:,1),y132(:,1)+(ones(1,length(y132))*hp)','m-');
legend('Ts = 10','Ts = 20','Ts = 25','Uklad ciagly');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Skok dodatni FD dla roznych ukladow dyskretnych - wysokosc');

%Skok TD
figure(25)
stairs(t254(:,1),y254(:,1)+(ones(1,length(y254))*Tp)');
hold on;
stairs(t264(:,1),y264(:,1)+(ones(1,length(y264))*Tp)');
hold on;
stairs(t274(:,1),y274(:,1)+(ones(1,length(y274))*Tp)');
hold on;
plot(t192(:,1),y192(:,1)+(ones(1,length(y192))*Tp)','m-');
legend('Ts = 10','Ts = 20','Ts = 25','Uklad ciagly');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Skok dodatni TD dla roznych ukladow dyskretnych - temperatura');

figure(26)
stairs(t254(:,1),y254(:,1)+(ones(1,length(y254))*hp)');
hold on;
stairs(t264(:,1),y264(:,1)+(ones(1,length(y264))*hp)');
hold on;
stairs(t274(:,1),y274(:,1)+(ones(1,length(y274))*hp)');
hold on;
plot(t192(:,1),y192(:,1)+(ones(1,length(y192))*hp)','m-');
legend('Ts = 10','Ts = 20','Ts = 25','Uklad ciagly');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Skok dodatni TD dla roznych ukladow dyskretnych - temperatura');

%%Wykresy bledow.
%Skok FC.
%Skok dodatni.
[t_count1,y_count1] = countErrorFunction(t11(1,:),y11(1,:),t12(:,1)',(y12(:,1)+(ones(1,length(y12))*Tp)')');
[t_count2,y_count2] = countErrorFunction(t21(1,:),y21(1,:),t22(:,1)',(y22(:,1)+(ones(1,length(y22))*Tp)')');
[t_count3,y_count3] = countErrorFunction(t31(1,:),y31(1,:),t32(:,1)',(y32(:,1)+(ones(1,length(y32))*Tp)')');
figure(27)
plot(t_count1,y_count1,'r-',t_count2,y_count2,'b-',t_count3,y_count3,'g-');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');legend('Skok o 5','Skok o 10','Skok o 15');
title('Wykres bledow linearyzacji przy skoku FH - temperatura');

[t_count4,y_count4] = countErrorFunction(t11(1,:),y11(2,:),t12(:,1)',(y12(:,2)+(ones(1,length(y12))*hp)')');
[t_count5,y_count5] = countErrorFunction(t21(1,:),y21(2,:),t22(:,1)',(y22(:,2)+(ones(1,length(y22))*hp)')');
[t_count6,y_count6] = countErrorFunction(t31(1,:),y31(2,:),t32(:,1)',(y32(:,2)+(ones(1,length(y32))*hp)')');
figure(28)
plot(t_count4,y_count4,'r-',t_count5,y_count5,'b-',t_count6,y_count6,'g-');
legend('Skok o 5','Skok o 10');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Wykres bledow linearyzacji przy skoku FH - wysokosc');

%Skok ujemny
[t_count7,y_count7] = countErrorFunction(t41(1,:),y41(1,:),t42(:,1)',(y42(:,1)+(ones(1,length(y42))*Tp)')');
[t_count8,y_count8] = countErrorFunction(t51(1,:),y51(1,:),t52(:,1)',(y52(:,1)+(ones(1,length(y52))*Tp)')');
[t_count9,y_count9] = countErrorFunction(t61(1,:),y61(1,:),t62(:,1)',(y62(:,1)+(ones(1,length(y62))*Tp)')');
figure(29)
plot(t_count7,y_count7,'r-',t_count8,y_count8,'b-',t_count9,y_count9,'g-');
legend('Skok o 5','Skok o 10','Skok o 15');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Wykres bledow linearyzacji przy skoku FH - temperatura');

[t_count10,y_count10] = countErrorFunction(t41(1,:),y41(2,:),t42(:,1)',(y42(:,2)+(ones(1,length(y42))*hp)')');
[t_count11,y_count11] = countErrorFunction(t51(1,:),y51(2,:),t52(:,1)',(y52(:,2)+(ones(1,length(y52))*hp)')');
[t_count12,y_count12] = countErrorFunction(t61(1,:),y61(2,:),t62(:,1)',(y62(:,2)+(ones(1,length(y62))*hp)')');
figure(28)
plot(t_count10,y_count10,'r-',t_count11,y_count11,'b-',t_count12,y_count12,'g-');
legend('Skok o 5','Skok o 10','Skok o 15');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Wykres bledow linearyzacji przy skoku FH - wysokosc');

%Skok FC.
%Skok dodatni.
[t_count1,y_count1] = countErrorFunction(t71(1,:),y71(1,:),t72(:,1)',(y72(:,1)+(ones(1,length(y72))*Tp)')');
[t_count2,y_count2] = countErrorFunction(t81(1,:),y81(1,:),t82(:,1)',(y82(:,1)+(ones(1,length(y82))*Tp)')');
[t_count3,y_count3] = countErrorFunction(t91(1,:),y91(1,:),t92(:,1)',(y92(:,1)+(ones(1,length(y92))*Tp)')');
figure(27)
plot(t_count1,y_count1,'r-',t_count2,y_count2,'b-',t_count3,y_count3,'g-');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
legend('Skok o 5','Skok o 10','Skok o 15');
title('Wykres bledow linearyzacji przy skoku FC - temperatura');

[t_count4,y_count4] = countErrorFunction(t71(1,:),y71(2,:),t72(:,1)',(y72(:,2)+(ones(1,length(y72))*hp)')');
[t_count5,y_count5] = countErrorFunction(t81(1,:),y81(2,:),t82(:,1)',(y82(:,2)+(ones(1,length(y82))*hp)')');
[t_count6,y_count6] = countErrorFunction(t91(1,:),y91(2,:),t92(:,1)',(y92(:,2)+(ones(1,length(y92))*hp)')');
figure(28)
plot(t_count4,y_count4,'r-',t_count5,y_count5,'b-',t_count6,y_count6,'g-');
legend('Skok o 5','Skok o 10','Skok o 15');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Wykres bledow linearyzacji przy skoku FC - wysokosc');

%Skok ujemny
[t_count7,y_count7] = countErrorFunction(t101(1,:),y101(1,:),t102(:,1)',(y102(:,1)+(ones(1,length(y102))*Tp)')');
[t_count8,y_count8] = countErrorFunction(t111(1,:),y111(1,:),t112(:,1)',(y112(:,1)+(ones(1,length(y112))*Tp)')');
[t_count9,y_count9] = countErrorFunction(t121(1,:),y121(1,:),t122(:,1)',(y122(:,1)+(ones(1,length(y122))*Tp)')');
figure(29)
plot(t_count7,y_count7,'r-',t_count8,y_count8,'b-',t_count9,y_count9,'g-');
legend('Skok o 5','Skok o 10','Skok o 15');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Wykres bledow linearyzacji przy skoku FC - temperatura');

[t_count10,y_count10] = countErrorFunction(t101(1,:),y101(2,:),t102(:,1)',(y102(:,2)+(ones(1,length(y102))*hp)')');
[t_count11,y_count11] = countErrorFunction(t111(1,:),y111(2,:),t112(:,1)',(y112(:,2)+(ones(1,length(y112))*hp)')');
[t_count12,y_count12] = countErrorFunction(t121(1,:),y121(2,:),t122(:,1)',(y122(:,2)+(ones(1,length(y122))*hp)')');
figure(30)
plot(t_count10,y_count10,'r-',t_count11,y_count11,'b-',t_count12,y_count12,'g-');
legend('Skok o 5','Skok o 10','Skok o 15');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Wykres bledow linearyzacji przy skoku FC - wysokosc');

%Skok FD.
%Skok dodatni.
[t_count1,y_count1] = countErrorFunction(t131(1,:),y131(1,:),t132(:,1)',(y132(:,1)+(ones(1,length(y132))*Tp)')');
[t_count2,y_count2] = countErrorFunction(t141(1,:),y141(1,:),t142(:,1)',(y142(:,1)+(ones(1,length(y142))*Tp)')');
[t_count3,y_count3] = countErrorFunction(t151(1,:),y151(1,:),t152(:,1)',(y152(:,1)+(ones(1,length(y152))*Tp)')');
figure(31)
plot(t_count1,y_count1,'r-',t_count2,y_count2,'b-',t_count3,y_count3,'g-');
legend('Skok o 5','Skok o 10','Skok o 15');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Wykres bledow linearyzacji przy skoku FD - temperatura');

[t_count4,y_count4] = countErrorFunction(t131(1,:),y131(2,:),t132(:,1)',(y132(:,2)+(ones(1,length(y132))*hp)')');
[t_count5,y_count5] = countErrorFunction(t141(1,:),y141(2,:),t142(:,1)',(y142(:,2)+(ones(1,length(y142))*hp)')');
[t_count6,y_count6] = countErrorFunction(t151(1,:),y151(2,:),t152(:,1)',(y152(:,2)+(ones(1,length(y152))*hp)')');
figure(32)
plot(t_count4,y_count4,'r-',t_count5,y_count5,'b-',t_count6,y_count6,'g-');
legend('Skok o 1','Skok o 2','Skok o 3');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Wykres bledow linearyzacji przy skoku FD - wysokosc');

%Skok ujemny
[t_count7,y_count7] = countErrorFunction(t161(1,:),y161(1,:),t162(:,1)',(y162(:,1)+(ones(1,length(y162))*Tp)')');
[t_count8,y_count8] = countErrorFunction(t171(1,:),y171(1,:),t172(:,1)',(y172(:,1)+(ones(1,length(y172))*Tp)')');
[t_count9,y_count9] = countErrorFunction(t181(1,:),y181(1,:),t182(:,1)',(y182(:,1)+(ones(1,length(y182))*Tp)')');
figure(33)
plot(t_count7,y_count7,'r-',t_count8,y_count8,'b-',t_count9,y_count9,'g-');
legend('Skok o 5','Skok o 10','Skok o 15');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Wykres bledow linearyzacji przy skoku FD - temperatura');

[t_count10,y_count10] = countErrorFunction(t161(1,:),y161(2,:),t162(:,1)',(y162(:,2)+(ones(1,length(y162))*hp)')');
[t_count11,y_count11] = countErrorFunction(t171(1,:),y171(2,:),t172(:,1)',(y172(:,2)+(ones(1,length(y172))*hp)')');
[t_count12,y_count12] = countErrorFunction(t181(1,:),y181(2,:),t182(:,1)',(y182(:,2)+(ones(1,length(y182))*hp)')');
figure(34)
plot(t_count10,y_count10,'r-',t_count11,y_count11,'b-',t_count12,y_count12,'g-');
legend('Skok o 5','Skok o 10','Skok o 15');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Wykres bledow linearyzacji przy skoku FD - wysokosc');

%Skok TD.
%Skok dodatni.
[t_count1,y_count1] = countErrorFunction(t191(1,:),y191(1,:),t192(:,1)',(y192(:,1)+(ones(1,length(y192))*Tp)')');
[t_count2,y_count2] = countErrorFunction(t201(1,:),y201(1,:),t202(:,1)',(y202(:,1)+(ones(1,length(y202))*Tp)')');
[t_count3,y_count3] = countErrorFunction(t211(1,:),y211(1,:),t212(:,1)',(y212(:,1)+(ones(1,length(y212))*Tp)')');
figure(35)
plot(t_count1,y_count1,'r-',t_count2,y_count2,'b-',t_count3,y_count3,'g-');
legend('Skok o 5','Skok o 10','Skok o 15');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Wykres bledow linearyzacji przy skoku TD - temperatura');

[t_count4,y_count4] = countErrorFunction(t191(1,:),y191(2,:),t192(:,1)',(y192(:,2)+(ones(1,length(y192))*hp)')');
[t_count5,y_count5] = countErrorFunction(t201(1,:),y201(2,:),t202(:,1)',(y202(:,2)+(ones(1,length(y202))*hp)')');
[t_count6,y_count6] = countErrorFunction(t211(1,:),y211(2,:),t212(:,1)',(y212(:,2)+(ones(1,length(y212))*hp)')');
figure(36)
plot(t_count4,y_count4,'r-',t_count5,y_count5,'b-',t_count6,y_count6,'g-');
legend('Skok o 5','Skok o 10','Skok o 15');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Wykres bledow linearyzacji przy skoku TD - wysokosc');

%Skok ujemny
[t_count7,y_count7] = countErrorFunction(t221(1,:),y221(1,:),t222(:,1)',(y222(:,1)+(ones(1,length(y222))*Tp)')');
[t_count8,y_count8] = countErrorFunction(t231(1,:),y231(1,:),t232(:,1)',(y232(:,1)+(ones(1,length(y232))*Tp)')');
[t_count9,y_count9] = countErrorFunction(t241(1,:),y241(1,:),t242(:,1)',(y242(:,1)+(ones(1,length(y242))*Tp)')');
figure(37)
plot(t_count7,y_count7,'r-',t_count8,y_count8,'b-',t_count9,y_count9,'g-');
legend('Skok o 5','Skok o 10','Skok o 15');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Wykres bledow linearyzacji przy skoku TD - temperatura');

[t_count10,y_count10] = countErrorFunction(t221(1,:),y221(2,:),t222(:,1)',(y222(:,2)+(ones(1,length(y222))*hp)')');
[t_count11,y_count11] = countErrorFunction(t231(1,:),y231(2,:),t232(:,1)',(y232(:,2)+(ones(1,length(y232))*hp)')');
[t_count12,y_count12] = countErrorFunction(t241(1,:),y241(2,:),t242(:,1)',(y242(:,2)+(ones(1,length(y242))*hp)')');
figure(38)
plot(t_count10,y_count10,'r-',t_count11,y_count11,'b-',t_count12,y_count12,'g-');
legend('Skok o 5','Skok o 10','Skok o 15');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Wykres bledow linearyzacji przy skoku TD - wysokosc');

%Porównanie transmitancji
%Porównanie transmitancji oraz modelu w przestrzeni stanow.
%Skok FH
figure(39);
plot(t12(:,1),y12(:,1)+(ones(1,length(y12))*Tp)','r-',t281(:,1),y281(:,1)+(ones(1,length(y281))*Tp)','b-',t291(:,1),y291(:,1)+(ones(1,length(y291))*Tp)','g-',t301(:,1),y301(:,1)+(ones(1,length(y301))*Tp)','m-', t311(:,1),y311(:,1)+(ones(1,length(y311))*Tp)','y-');
legend('równania stanu','transmitancja ciagla','transmitancja dyskretna Ts = 10', 'transmitancja dyskretna Ts = 20', 'transmitancja dyskretna Ts = 25');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Porównanie modeli dyskretnych przy skoku FH - temperatura');

figure(40);
plot(t12(:,1),y12(:,2)+(ones(1,length(y12))*hp)','r-',t281(:,1),y281(:,2)+(ones(1,length(y281))*hp)','b-',t291(:,1),y291(:,2)+(ones(1,length(y291))*hp)','g-',t301(:,1),y301(:,2)+(ones(1,length(y301))*hp)','m-', t311(:,1),y311(:,2)+(ones(1,length(y311))*hp)','y-');
legend('równania stanu','transmitancja ciagla','transmitancja dyskretna Ts = 10', 'transmitancja dyskretna Ts = 20', 'transmitancja dyskretna Ts = 25');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Porównanie modeli dyskretnych przy skoku FH - wysokosc');

%Skok FC
figure(41);
plot(t72(:,1),y72(:,1)+(ones(1,length(y72))*Tp)','r-',t282(:,1),y282(:,1)+(ones(1,length(y282))*Tp)','b-',t292(:,1),y292(:,1)+(ones(1,length(y292))*Tp)','g-',t302(:,1),y302(:,1)+(ones(1,length(y302))*Tp)','m-', t312(:,1),y312(:,1)+(ones(1,length(y312))*Tp)','y-');
legend('równania stanu','transmitancja ciagla','transmitancja dyskretna Ts = 10', 'transmitancja dyskretna Ts = 20', 'transmitancja dyskretna Ts = 25');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Porównanie modeli dyskretnych przy skoku FC - temperatura');

figure(42);
plot(t72(:,1),y72(:,2)+(ones(1,length(y72))*hp)','r-',t282(:,1),y282(:,2)+(ones(1,length(y282))*hp)','b-',t292(:,1),y292(:,2)+(ones(1,length(y292))*hp)','g-',t302(:,1),y302(:,2)+(ones(1,length(y302))*hp)','m-', t312(:,1),y312(:,2)+(ones(1,length(y312))*hp)','y-');
legend('równania stanu','transmitancja ciagla','transmitancja dyskretna Ts = 10', 'transmitancja dyskretna Ts = 20', 'transmitancja dyskretna Ts = 25');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Porównanie modeli dyskretnych przy skoku FC - wysokosc');

%Skok FD
figure(43);
plot(t132(:,1),y132(:,1)+(ones(1,length(y132))*Tp)','r-',t283(:,1),y283(:,1)+(ones(1,length(y283))*Tp)','b-',t293(:,1),y293(:,1)+(ones(1,length(y293))*Tp)','g-',t303(:,1),y303(:,1)+(ones(1,length(y303))*Tp)','m-', t313(:,1),y313(:,1)+(ones(1,length(y313))*Tp)','y-');
legend('równania stanu','transmitancja ciagla','transmitancja dyskretna Ts = 10', 'transmitancja dyskretna Ts = 20', 'transmitancja dyskretna Ts = 25');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Porównanie modeli dyskretnych przy skoku FD - temperatura');

figure(44);
plot(t132(:,1),y132(:,2)+(ones(1,length(y132))*hp)','r-',t283(:,1),y283(:,2)+(ones(1,length(y283))*hp)','b-',t293(:,1),y293(:,2)+(ones(1,length(y293))*hp)','g-',t303(:,1),y303(:,2)+(ones(1,length(y303))*hp)','m-', t313(:,1),y313(:,2)+(ones(1,length(y313))*hp)','y-');
legend('równania stanu','transmitancja ciagla','transmitancja dyskretna Ts = 10', 'transmitancja dyskretna Ts = 20', 'transmitancja dyskretna Ts = 25');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Porównanie modeli dyskretnych przy skoku FD - wysokosc');

%Skok TD
figure(45);
plot(t192(:,1),y192(:,1)+(ones(1,length(y192))*Tp)','r-',t284(:,1),y284(:,1)+(ones(1,length(y284))*Tp)','b-',t294(:,1),y294(:,1)+(ones(1,length(y294))*Tp)','g-',t304(:,1),y304(:,1)+(ones(1,length(y304))*Tp)','m-', t314(:,1),y314(:,1)+(ones(1,length(y314))*Tp)','y-');
legend('równania stanu','transmitancja ciagla','transmitancja dyskretna Ts = 10', 'transmitancja dyskretna Ts = 20', 'transmitancja dyskretna Ts = 25');
xlabel('czas [s]');
ylabel('Temperatura [st. C]');
title('Porównanie modeli dyskretnych przy skoku TD - temperatura');

figure(46);
plot(t192(:,1),y192(:,2)+(ones(1,length(y192))*hp)','r-',t284(:,1),y284(:,2)+(ones(1,length(y284))*hp)','b-',t294(:,1),y294(:,2)+(ones(1,length(y294))*hp)','g-',t304(:,1),y304(:,2)+(ones(1,length(y304))*hp)','m-', t314(:,1),y314(:,2)+(ones(1,length(y314))*hp)','y-');
legend('równania stanu','transmitancja ciagla','transmitancja dyskretna Ts = 10', 'transmitancja dyskretna Ts = 20', 'transmitancja dyskretna Ts = 25');
xlabel('czas [s]');
ylabel('Wysokosc [cm]');
title('Porównanie modeli dyskretnych przy skoku TD - wysokosc');

