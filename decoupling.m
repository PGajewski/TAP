clear all;
close all;

%wpisanie danych
C=0.4;
a=8;
global THp;
global TCp;
TCp=21;
THp=64;
TDp=30;
FCp=32;
FHp=20;
FDp=9;
tau=40;
tau_c=100;
hp=58.14;
Tp=36.43;
Tsim = 1000;

%Obliczanie transmitancji
%Wyznaczanie wartoœci liczbowych macierzy A
a11=-(FHp+FCp+FDp)/(C*hp^2);
a12=(-2)*((THp-Tp)*FHp+(TCp-Tp)*FCp+(TDp-Tp)*FDp)/(C*hp^3);
a21=0;
a22=((FHp+FCp+FDp)/(hp^2)-a/(2*(sqrt(hp^3))))/(-2*C);

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
c=[1 0; 0 1];

%Wyznaczanie wartoœci macierzy D
D=[0 0 0 0; 0 0 0 0];

%Obliczenie transmitancji
syms s;
tf('s');
G1=c*(inv(s*eye(2)-A))*B+D;
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

%%
% G1=c2d(G1,1);
D12=-G1(1,2)/G1(1,1);
[num,den]=tfdata(D12);  %dodanie opoznienia do transmitancji GhFC
T=1;

D12 = tf(num,den,'InputDelay',tau_c);
D12_z = c2d(D12, T);

D21=-G1(2,1)/G1(2,2);
D21_z = c2d(D21, T);
%dodawanie opóŸnieñ do transmitancji
%G111 = tf(G1(1,2),'InputDelay',tauC)

%Create PIDs.
Kp1 = 0.600392409469819;
Ki1 = 0.0259650078715821;
Kd1 = 0;
N1 = 100;

Kp2 = 0.180678962603669;
Ki2 = 0.00106729417698503;
Kd2 = 0;
N2 = 100;

pid1 = pid(Kp1, Ki1, Kd1, N1, T);
pid2 = pid(Kp2, Ki2, Kd2, N2, T);

Tzad_change = 100;
hzad_change = 500;

%Set values.
Tzad = [ones(1,Tzad_change/T)*Tp, ones(1,(Tsim-Tzad_change)/T)*(Tp+3)];
hzad = [ones(1,hzad_change/T)*Tp, ones(1,(Tsim-hzad_change)/T)*(hp+3)];

FH = [];
FC = [];
x_sim = [];
t_sim = [];
Tout = zeros(1,length(1:1:Tsim/T));
hout = zeros(1,length(1:1:Tsim/T));
actual_T = Tp;
actual_h = hp;
pid1_x = [FHp];
pid2_x = [FCp];
decoupler12_x = [0];
decoupler21_x = [0];
x_0 = [Tp, hp];
for i=1:1:Tsim/T
    i
    %Count control signals.
    e1 = Tzad(i) - actual_T;
    e2 = hzad(i) - actual_h;
    
    %Count PID outputs.
    [y1, t, pid1_x] = lsim(pid1, [e1; e1], [0; T], pid1_x);
    [y2, t, pid2_x] = lsim(pid2, [e2; e2], [0; T], pid2_x);
    
    %Count decouplers outputs.
    [d1, t, decoupler12_x] = lsim(D12_z, y2, [0; T], decoupler12_x);
    [d2, t, decoupler21_x] = lsim(D21_z, y1, [0; T], decoupler21_x);
    
    %Count input of object.
    u1 = y1(2,1) + d1(2,1);
    u2 = y2(2,1) + d2(2,1);
    
    %Choose FC and FH
    FH = [FH, u1];
    FC = [FC, u2];
    if i*T <= tau_c
        Fcin = FCp;
    else
        Fcin = FC(i-tau_c/T);
    end
    Fhin = u1;

    %Simulate object.
    stateHandler = @(t,x) stateFunction(t,x,Fhin, Fcin, FDp, THp, TCp, TDp);
    [temp1,temp2]=ode45(stateHandler,[0 T],x_0);
    t_sim = [t_sim temp1(end,end)];
    x_sim = [x_sim temp2(end,:)'];
    x_0 = temp2(end,:);
    
    %Choose outputs.
    if i*T <= tau
        actual_T = Tp;
    else
        actual_T = x_sim(1, i-tau/T);
    end
    actual_h = x_sim(2, end);
    
    %Save outputs.
    Tout(i) = actual_T;
    hout(i) = actual_h;
end
t = (T:T:Tsim);

figure(1);
plot(t,Tout,t,hout);

    