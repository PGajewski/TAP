function [F] = flowOutput(h)
%flowOutput Toriccelli equation of flow from the tank.
alfa = 8;
F = alfa * sqrt(h);
end

