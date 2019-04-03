function [x_next] = stateFunction(t,x,F_H, F_C, F_D, T_H, T_C, T_D)
%stateFunction - equation of state in model.
C = 0.4;
alfa = 8;
h_next = (1/(2*C*x(2)))*(F_H + F_C + F_D - alfa*sqrt(x(2)));
T_next = (1/(C*(x(2)^2)))*((T_H-x(1))*F_H + (T_C - x(1))*F_C + (T_D-x(1))*F_D);
x_next = [T_next;h_next];
end

