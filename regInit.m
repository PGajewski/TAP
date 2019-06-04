%classPID(K, Ti, Kd, Td, Tp, Hlim, Llim, Dir, ManVal)

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

Ts = 1;
Kp1 = 0.600392409469819;
Ti1 = 38.51337172496948;
%Ti1 = 100;
Kd1 = 0;
Td1 = 0;

Kp2 = 0.180678962603669;
%Ki2 = 936.950914688;
Ti2 = 50;
Td2 = 0;
Td2 = 0;

global reg1;
global reg2;
reg1=classPID(Kp1, Ti1, Kd1, Td1, Ts, 50, 0, 1, 1);
reg2=classPID(Kp2, Ti2, Kd2, Td2, Ts, 50, 0, 1, 1);

%SetAutoMan(obj, AutoMan, ManVal)

reg1.SetAutoMan(1,FHp);
reg1.stateUpdate(FHp);

reg2.SetAutoMan(1,FCp);
reg2.stateUpdate(FCp);
