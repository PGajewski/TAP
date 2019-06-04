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
Ts = 1;
Tsim = 100;

%Badane K_kryt
K_kryt = 1;

%D regulator
%classPID(K, Ti, Kd, Td, Tp, Hlim, Llim, Dir, ManVal) 
reg=classPID(0, 0, 1, 30, Ts, 50, 0, 1, 1);

%SetAutoMan(obj, AutoMan, ManVal)
reg.SetAutoMan(1,1);

%% Regulator FH/T
%PID
%reTune(obj, K, Ti, Kd, Td)
reg.reTune(K_kryt, 9999, 0, 0);

stpt = FHp;
pv=Tp;
u=0;

t_sim = [];
x_sim = [];
y_sim = [];

x_0 = [Tp, hp];
for i=1:Ts:Tsim
  u = reg.calc(pv,stpt)
  
  %Prepare ODE options.
  options = odeset('RelTol',1e-8,'AbsTol',1e-10);
  
  %wyliczenie wyjscia obiektu na nastepny krok wykorzystujac u
  stateHandler = @(t,x) stateFunction(t, x, u , FCp, FDp, THp, TCp, TDp);
  [temp1,temp2]=ode45(stateHandler,[(i-1) * Ts i * Ts],x_0, options);
  
  t_sim = [t_sim temp1(1:end-1,:)'];
  x_sim = [x_sim temp2(1:end-1,:)'];
  x_0 = x_sim(:, end);
  
    %Wyliczenie wyjscia.
  t = (i-1) * Ts;
  if t < tau
      pv = Tp;
  else
      %Count expected process time for get temperature.
      expected_time = t - tau;
      index = (expected_time+1)/Ts;
      pv = x_sim(1,index);
  end
  
  y_sim = [y_sim, pv];
end

figure(1);
plot(0:Ts:Tsim-1,y_sim(1,:));
legend('T');
xlabel('Czas [s]');
ylabel('Temperatura');
title('Przebieg temperatury (symulacja)');

%% Regulator FC/H
%PID
%reTune(obj, K, Ti, Kd, Td)
reg.reTune(K_kryt, 9999, 0, 0);

stpt = FCp;
pv=hp;
u=0;

t_sim = [];
x_sim = [];
y_sim = [];
ctr_history = ones(1,length(1:Ts:Tsim));
x_0 = [Tp, hp];
for i=1:Ts:Tsim
  u = reg.calc(pv,stpt)
  ctr_history(1,i) = u;
      %Wyliczenie wyjscia.
  t = (i-1) * Ts;
  if t < tau_C
      u = TCp;
  else
      %Count expected process time for get temperature.
      expected_time = t - tau_C;
      index = (expected_time+1)/Ts;
      u = ctr_history(1,index);
  end
  %Prepare ODE options.
  options = odeset('RelTol',1e-8,'AbsTol',1e-10);
  
  %wyliczenie wyjscia obiektu na nastepny krok wykorzystujac u
  stateHandler = @(t,x) stateFunction(t, x, FHp , u, FDp, THp, TCp, TDp);
  [temp1,temp2]=ode45(stateHandler,[(i-1) * Ts i * Ts],x_0, options);
  
  t_sim = [t_sim temp1(1:end-1,:)'];
  x_sim = [x_sim temp2(1:end-1,:)'];
  x_0 = x_sim(:, end);
  
  %Wyliczenie wyjœcia
  pv = x_sim(2,end);
  
  y_sim = [y_sim, pv];
end

figure(2);
plot(0:Ts:Tsim-1,y_sim(1,:));
legend('H');
xlabel('Czas [s]');
ylabel('Wysokosc');
title('Przebieg wysokosci (symulacja)');
