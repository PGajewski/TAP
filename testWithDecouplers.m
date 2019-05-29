clear all;

%% Przygotowanie transmitancji
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

%dla D12
[numD12_z,denD12_z] = tfdata(D12_z);
numD12_z=numD12_z{1};
denD12_z=denD12_z{1};

%dla D21
[numD21_z, denD21_z]=tfdata(D21_z);
numD21_z=numD21_z{1};
denD21_z=denD21_z{1};


%% Regulacja
%D regulator
%classPID(K, Ti, Kd, Td, Tp, Hlim, Llim, Dir, ManVal) 
%D regulator
%p=classPID(0, 0, 1, 30, 1, 100, -100, 1, 1)
reg1=classPID(0, 0, 1, 30, 1, 100, -100, 1, 1);
reg2=classPID(0, 0, 1, 30, 1, 100, -100, 1, 1);

%SetAutoMan(obj, AutoMan, ManVal)
reg1.SetAutoMan(1,1)
reg2.SetAutoMan(1,1)

hp=58.140625;
Tp=36.4262295081967;

Ts = 1;
Kp1 = 1.600392409469819;
Ti1 = 45.51337172496948;
%Ti1 = 100;
Kd1 = 0;
Td1 = 0;

Kp2 = 0.180678962603669;
%Ki2 = 936.950914688;
%Ti2 = 5000;
Ti2=936.94880152431125234071143952989;
Td2 = 0;
Kd2 = 0;

%PID
%reTune(obj, K, Ti, Kd, Td)
reg1.reTune(Kp1, Ti1, Kd1, Td1);
reg2.reTune(Kp2, Ti2, Kd2, Td2);

% stpt = 1;
% pv=0;
% u=0;
FH=ones(5000,1)*FHp;
FCin=ones(5000,1)*FCp;
h=ones(5000,1)*hp;
T=ones(5000,1)*Tp;
Tout=ones(5000,1)*Tp;

d12=zeros(1,5000);
d21=zeros(1,5000);

TpSetPoint=[ones(1,4000)* Tp ones(1,6000)*(Tp+1)];

%obiekt = tf(1,[30 1]);
loopSize=10000;
for i=200:1:loopSize
% 
  FH(i) = reg1.calc(Tout(i-1),TpSetPoint(i));
% 
  FCin(i) = reg2.calc(h(i-1),hp);
  
  if i < 200
        FH(i)=FHp;
        FCin(i)=FCp;    
  end
  
    %Count decouplers outputs.
    %dla D12
    d12(i)=-denD12_z(2)*d12(i-1)-denD12_z(3)*d12(i-2)-denD12_z(4)*d12(i-3)+numD12_z(1)*FH(i-tau_c)+numD12_z(2)*FH(i-1-tau_c)+numD12_z(3)*FH(i-2-tau_c)+numD12_z(4)*FH(i-3-tau_c);
    d21(i)=numD21_z*FCin(i);
    FHd12=FH(i)+d21(i);
    FCd21=FCin(i)+d12(i);

    FH(i)=FHd12 + FHp;
    FCin(i)=FCd21 + FCp;

        
  %wyliczenie wyjscia obiektu na nastepny krok wykorzystujac u
  FC(i)=FCin(i-tau_c);
  h(i)=h(i-1)+1/(2*C*h(i-1))*(FH(i)+FC(i)+FDp-a*sqrt(h(i-1)));
  T(i)=T(i-1)+1/(C*(h(i-1))^2)*((THp-T(i-1))*FH(i)+(TCp-T(i-1))*FC(i)+(TDp-T(i-1))*FDp);
  Tout(i)=T(i-tau);
end

figure();
plot(h(200:loopSize));
figure();
plot(Tout(200:loopSize));

figure();
plot(d12(200:loopSize));

