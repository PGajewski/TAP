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
Ts =1;
Tsim = 2000;

%Obliczanie transmitancji
% %Wyznaczanie wartosci liczbowych macierzy A
% a11=-(FHp+FCp+FDp)/(C*hp^2);
% a12=(-2)*((THp-Tp)*FHp+(TCp-Tp)*FCp+(TDp-Tp)*FDp)/(C*hp^3);
% a21=0;
% a22=((FHp+FCp+FDp)/(hp^2)-a/(2*(sqrt(hp^3))))/(-2*C);
% sys
% A=[a11 a12; a21 a22];
% 
% %Wyznaczanie wartosci liczbowych macierzy B
% b11=1/(C*hp^2)*(THp-Tp);
% b12=1/(C*hp^2)*(TCp-Tp);
% b13=1/(C*hp^2)*(TDp-Tp);
% b14=1/(C*hp^2)*FDp;
% b21=1/(2*C*hp);
% b22=1/(2*C*hp);
% b23=1/(2*C*hp);
% b24=0;
% 
% B=[b11 b12 b13 b14; b21 b22 b23 b24];
% 
% %Wyznaczanie wartosci liczbowych macierzy C
% C=[1 0; 0 1];
% 
% %Wyznaczanie wartosci macierzy D
% D=[0 0 0 0; 0 0 0 0];
% 
% sys_con = ss(A,B,C,D,'InputDelay',[0,tau_C,0,0],'OutputDelay',[tau,0]);
trans = tf([1], [30 1]);
sys_con = ss(trans);
sys_dyskr = c2d(sys_con, Ts);

N = 100;
Nu = 10;

n = size(sys_dyskr.A);
teta = obsv(sys_dyskr);

if rank(teta) == n(1,1)
    
else
    Kobsp = place(sys_dyskr.A', sys_dyskr.C', zeros(1, n(1,1)));
    Kobsb = A^-1 * Kobsp;
end

Y_zad = [];
u = [];
[M,CtAt,CtV]=MPCSmatrices(sys_dyskr.A, sys_dyskr.B, sys_dyskr.C, N, Nu);

for k = 1 : Tsim
    delta_u = K*(Y_zad(k) - CtAt - CtV * sys.dyskr.B * u(k-1) - CtV*v(k);
    u(k) = u(k-1) + delta_u;
end