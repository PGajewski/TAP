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
Tsim = 2000;
lambda_coefficient = 1;

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
Nu = 50;

% restrictions on outputs.
xMin=[-inf;-inf];
xMax=[inf;inf];



teta = obsv(sys_dyskr);
A_dyskr =sys_dyskr.A;
B1_dyskr= sys_dyskr.B(:,1:2);
B2_dyskr= sys_dyskr.B(:,3:4);
C_dyskr=sys_dyskr.C;
D_dyskr=sys_dyskr.D;
[M,CtAt,CtV]=MPCSmatrices(A_dyskr, B1_dyskr, C_dyskr, N, Nu);

[nx, nu]=size(B1_dyskr); ny=size(C_dyskr,1);

psi=eye(N*ny);
lambda=eye(Nu*nu)*lambda_coefficient;

% restrictions on control signal.
uMin=ones(Nu*nu,1) * (-5);
deltaUMin=ones(nu*Nu,1) * (-0.5);
deltaUMax=ones(nu*Nu,1) * (0.5);
uMax=ones(Nu*nu,1) * 5;
y_min_val= ones(N*ny,1) * (-10); %*. [-Tp;-hp];
y_max_val= ones(N*ny,1) * 10 ;%*. [100-Tp;100-hp];

%%Simulation
sim_len = length(1:Ts:Tsim);
u=zeros(nu,sim_len);
U=zeros(Nu*nu,sim_len);
y=zeros(ny,sim_len);
x=zeros(nx,sim_len);
v=zeros(nx,sim_len);
y_zad=[zeros(ny,1000),ones(ny,sim_len-1000)];
%y_zad=[zeros(ny,1000),zeros(ny,sim_len-1000)];
FD = ones(1,sim_len)*FDp;
TD = ones(1,sim_len)*TDp;

h = ones(1,sim_len)*hp;
T = ones(1,sim_len)*Tp;
Tout = ones(1,sim_len)*Tp;

%"Offline" matrices.
H_quad = 2 * (M'*psi*M + lambda);
J = [];
for i=1:Nu
    W = [];
    for j=1:i
       W = [W eye(nu)]; 
    end
    W = [W zeros(nu,nu * (Nu-i))];
    J = [J;W];
end
A_quad = [-J;J;-M;M];

options_quad = optimoptions('quadprog','Display','off','MaxIterations', 400)
delta_U = zeros(1,Nu);
for k = 200 : sim_len-N
    k
    %Prepare y_zad.
    y_zad_act = zeros(1,N*ny);
    for i=2:2:2*N
       y_zad_act(1,i-1) = y_zad(1,k-1+i/2);
       y_zad_act(1,i) = y_zad(2,k-1+i/2);
    end
    
    %Error.
    v(:,k)=x(:,k)-(A_dyskr*x(:,k-1)+B1_dyskr*[u(1,k-1);u(2,k-1-floor(k-tau_C/Ts))] + B2_dyskr*[FD(1,k-1)-FDp;TD(1,k-1)-TDp]);
    y0 = CtAt*[x(1,k-floor(tau/Ts)); x(2,k)] + CtV * [B1_dyskr, B2_dyskr] * [u(1,k-1);u(2,k-1-floor(k-tau_C/Ts));FD(1,k)-FDp;TD(1,k)-TDp] + CtV*v(:,k);

    f_quad = -2 * M' * psi * (y_zad_act' - y0);
    
    u_min_res = -uMin + U(:,k-1);
    u_max_res = uMax - U(:,k-1);
    y_min_res = -y_min_val + y0;
    y_max_res = y_max_val -y0;
    b_quad = [u_min_res; u_max_res; y_min_res; y_max_res];
    
    [delta_U] = quadprog(H_quad,f_quad,A_quad,b_quad, [], [], deltaUMin, deltaUMax, delta_U, options_quad);
    %[delta_U] = quadprog(H_quad,f_quad,A_quad,b_quad, [], [], [], [], delta_U, options_quad);
    U(:,k) = U(:,k-1)+delta_U;

    u(:,k) = U(1:nu,k);
    
    %Count object.
  h(k+1)=h(k)+1/(2*C_coefficient*h(k))*(u(1,k)+ FHp +u(2,k-floor(k-tau_C/Ts)) +FCp +FD(1,k)-a*sqrt(h(k)));
  T(k+1)=T(k)+1/(C_coefficient*(h(k))^2)*((THp-T(k))*(u(1,k) + FHp)+(TCp-T(k))*(u(2,k-floor(k-tau_C/Ts)) + FCp)+(TD(1,k)-T(i))*FD(1,k));

    %Restrictions
    h(k+1)=min(h(k+1),xMax(1));
    h(k+1)=max(h(k+1),xMin(1));        
 
    T(k+1)=min(T(k+1),xMax(2));
    T(k+1)=max(T(k+1),xMin(2)); 
    
  Tout(k)=T(k-floor(tau/Ts));
  x(:,k+1) = [T(k+1)-Tp;h(k+1)-hp];
  y(:,k) = [Tout(k)-Tp;h(k)-hp];
    
end

figure(1);
stairs(1:sim_len-N, Tout(1:sim_len-N));

figure(2);
stairs(1:sim_len-N, h(1:sim_len-N));


