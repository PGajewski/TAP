clear all;
close all;

%wpisanie danych
C_coefficient=0.4;
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
Tsim = 3000;
lambda_coefficient = 0.001;
psi_coefficient = 1;

%Obliczanie transmitancji
%Wyznaczanie wartosci liczbowych macierzy A
a11=-(FHp+FCp+FDp)/(C_coefficient*hp^2);
a12=(-2)*((THp-Tp)*FHp+(TCp-Tp)*FCp+(TDp-Tp)*FDp)/(C_coefficient*hp^3);
a21=0;
a22=((FHp+FCp+FDp)/(hp^2)-a/(2*(sqrt(hp^3))))/(-2*C_coefficient);
A=[a11 a12; a21 a22];

%Wyznaczanie wartosci liczbowych macierzy B
b11=1/(C_coefficient*hp^2)*(THp-Tp);
b12=1/(C_coefficient*hp^2)*(TCp-Tp);
b13=1/(C_coefficient*hp^2)*(TDp-Tp);
b14=1/(C_coefficient*hp^2)*FDp;
b21=1/(2*C_coefficient*hp);
b22=1/(2*C_coefficient*hp);
b23=1/(2*C_coefficient*hp);
b24=0;

B=[b11 b12 b13 b14; b21 b22 b23 b24];

%Wyznaczanie wartosci liczbowych macierzy C
C=[1 0; 0 1];

%Wyznaczanie wartosci macierzy D
D=[0 0 0 0; 0 0 0 0];

sys_con = ss(A,B,C,D,'InputDelay',[0,tau_C,0,0],'OutputDelay',[tau,0]);
% trans = tf([1], [30 1]);
% sys_con = ss(trans);
sys_dyskr = c2d(sys_con, Ts);

N = 300;
Nu = 20;

% restrictions on control signal.
uMin=[-inf,-inf];
deltaUMin=[-inf,-inf];
deltaUMax=[inf,inf];
uMax=[inf, inf];

% restrictions on outputs.
xMin=[-inf,-inf];
xMax=[inf,inf];

teta = obsv(sys_dyskr);
A_dyskr =sys_dyskr.A;
B1_dyskr= sys_dyskr.B(:,1:2);
B2_dyskr= sys_dyskr.B(:,3:4);
C_dyskr=sys_dyskr.C;
D_dyskr=sys_dyskr.D;
[M,CtAt,CtV]=MPCSmatrices(A_dyskr, B1_dyskr, C_dyskr, N, Nu);

[nx, nu]=size(B1_dyskr); ny=size(C_dyskr,1);

psi=eye(N*ny)*psi_coefficient;
lambda=eye(Nu*nu)*lambda_coefficient;
K=(M'*psi*M+lambda)^-1*M'*psi;
K1=K(1:nu,:);

% if rank(teta) == nx
%     
% else
%     Kobsp = place(A_dyskr', C_dyskr', zeros(1, nx));
%     Kobsb = A^-1 * Kobsp;
% end

%%Simulation
sim_len = length(1:Ts:Tsim);
u=zeros(nu,sim_len-1000);
y=zeros(ny,sim_len);
x=zeros(nx,sim_len);
v=zeros(nx,sim_len);
%y_zad=[ones(ny,1000),ones(ny,sim_len-1000)];
%y_zad=[zeros(ny,1000),[ones(1,sim_len-1000)*0.5; zeros(1,sim_len-1000)]];

y_zad=[zeros(ny,1000),zeros(ny,sim_len-1000)];
FD = ones(1,sim_len)*FDp;
TD = ones(1,sim_len)*TDp;

h = ones(1,sim_len)*hp;
T = ones(1,sim_len)*Tp;
Tout = ones(1,sim_len)*Tp;


for k = 200 : sim_len-N
    %Prepare y_zad.
    y_zad_act = zeros(1,N*ny);
    for i=2:2:2*N
       y_zad_act(1,i-1) = y_zad(1,k-1+i/2);
       y_zad_act(1,i) = y_zad(2,k-1+i/2);
    end
    
%     %Prepare y_zad.
%     y_zad_act = zeros(1,N*ny);
%     for i=2:2:2*N
%        y_zad_act(1,i-1) = y_zad(1,k);
%        y_zad_act(1,i) = y_zad(2,k);
%     end
%     
    %Error.
%     v(:,k)=x(:,k)-(A_dyskr*x(:,k-1)+B1_dyskr*[u(1,k-1);u(2,k-1-floor(tau_C/Ts))] + B2_dyskr*[FD(1,k-1)-FDp;TD(1,k-1)-TDp]);
%     delta_u = K1*(0 - CtAt*x(:,k) - CtV * [B1_dyskr, B2_dyskr] * [u(1,k-1);u(2,k-1-floor(tau_C/Ts));FD(1,k)-FDp;TD(1,k)-TDp] - CtV*v(:,k));

    v(:,k)=x(:,k)-(A_dyskr*x(:,k-1)+B1_dyskr*[u(1,k-1);u(2,k-1-floor(tau_C/Ts))] + B2_dyskr*[FD(1,k-1)-FDp;TD(1,k-1)-TDp]);
    delta_u = K1*(y_zad_act' - CtAt*x(:,k) - CtV * B1_dyskr * [u(1,k-1);u(2,k-1-floor(tau_C/Ts))]);

    
    %Restrictions.
    for i=1:nu
        delta_u(i,1) = min(delta_u(i,1), deltaUMax(i));
        delta_u(i,1) = max(delta_u(i,1), deltaUMin(i));
    end
    delta_u
    u(:,k) = u(:,k-1) + delta_u;
    
    for i=1:nu
        u(i,k) = min(u(i,k), uMax(i));
        u(i,k) = max(u(i,k), uMin(i));
    end
    
    %Count object
    
  h(k+1)=h(k)+1/(2*C_coefficient*h(k))*(u(1,k)+ FHp +u(2,k-floor(tau_C/Ts)) +FCp +FD(1,k)-a*sqrt(h(k)));
  T(k+1)=T(k)+1/(C_coefficient*(h(k))^2)*((THp-T(k))*(u(1,k) + FHp)+(TCp-T(k))*(u(2,k-floor(tau_C/Ts)) + FCp)+(TD(1,k)-T(i))*FD(1,k));

    %Restrictions
%     h(k+1)=min(h(k+1),xMax(1));
%     h(k+1)=max(h(k+1),xMin(1));        
%  
%     T(k+1)=min(T(k+1),xMax(2));
%     T(k+1)=max(T(k+1),xMin(2)); 
    
  Tout(k)=T(k-floor(tau/Ts));
  x(:,k+1) = [T(k+1)-Tp;h(k+1)-hp];
    
end

figure(1);
stairs(1:sim_len-N, Tout(1:sim_len-N));
hold on;
stairs(1:sim_len-N, y_zad(2,1:sim_len-N)+Tp);
title('Analitycal MPCS - temperature');
xlabel('k')
ylabel('temperature T [C degrees]');

figure(2);
stairs(1:sim_len-N, h(1:sim_len-N));
hold on;
stairs(1:sim_len-N, y_zad(1,1:sim_len-N)+hp);
title('Analitycal MPCS - height');
xlabel('k')
ylabel('height [cm]');


