T_C=21; T_H=64; T_D=30;
F_C = 32; F_H = 20; F_D = 9;

stateHandler = @(t,x) stateFunction(t,x,F_H, F_C, F_D, T_H, T_C, T_D);
ode45(stateHandler,[0 10],[58.14 36.43]);