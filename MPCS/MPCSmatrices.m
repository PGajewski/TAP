function [M,CtAt,CtV]=MPCSmatrices(A,B,C,N,Nu)
%function "MPCSmatrices": macierze M, CtAt i CtV ("t"- tilda) regulatora MPCS,
%DANE WEJŒCIOWE (wymagane):
%macierze A, B i C równañ stanu obiektu: x(k+1)=Ax(k)+Bu(k), y(k)=Cx(k),
%N - horyzont predykcji, Nu - horyzont sterowania.
[nx nu]=size(B); ny=size(C,1);
%macierze CtAt, CtV (t-tilda) i M1 (pierwsza kolumna (blokowa) macierzy M):
X=eye(nx); Xa=A;
M1=C*B; CtAt=C*A; CtV=C;
for i=2:N
Xa=A*Xa;
CtAt=[CtAt;C*Xa];
X=eye(nx)+A*X;
M1=[M1; C*X*B];
CtV=[CtV;C*X];
end
%macierz dynamiczna M (wykorzystanie struktury macierzy Toeplitza):
M=M1;
for i=2:Nu
M=[M [zeros((i-1)*ny,nu); M1(1:(N-i+1)*ny,:)]];
end