function [t,x,y] = linearSimulation(FH_with_time,FC_with_time,FD_with_time, TD_with_time, G, Ts, sim_time, x_0)
%objectSimulation Simulation of object
%   FH_with_time - matrix of changing FH values: new value (first column) and time step (second column).
%   FC_with_time - matrix of changing FC values: new value (first column) and time step (second column).
%   FD_with_time - matrix of changing FD values: new value (first column) and time step (second column).
%   tau_H        - transport delay of signal FH.
%   tau_C        - transport delay of signal FC.
%   tau_D        - transport delay of signal FD.

% Check variables dimension.
temp1 = size(FH_with_time);
temp2 = size(FC_with_time);
temp3 = size(FD_with_time);
temp4 = size(TD_with_time);
temp5 = size(G);
assert((temp1(2) == 2 && temp2(2) == 2 && temp3(2) == 2 && temp4(2) == 2),'Incompatible changing moment values!');
assert((isvector(x_0) && length(x_0) == 2),'Init condition must be a vector with size 2'); 
assert((temp5(1) == 2 && temp5(2) == 4),'Bad transfer function dimension');

%Check last value change moment timestep.
assert((FH_with_time(end,1)<=sim_time) && (FC_with_time(end,1)<=sim_time) && (FD_with_time(end,1)<=sim_time) && (TD_with_time(end,1)<=sim_time),'One or more value changes is outside simulation range!');

FH_length = temp1(1)+1;
FC_length = temp2(1)+1;
FD_length = temp3(1)+1;
TD_length = temp4(1)+1;
counter_FC = 2;
counter_FH = 2;
counter_FD = 2;
counter_TD = 2;

%Add simulation end to changes moments.
FH = [FH_with_time;[sim_time,FH_with_time(end,2)]];
FC = [FC_with_time;[sim_time,FC_with_time(end,2)]];
FD = [FD_with_time;[sim_time,FD_with_time(end,2)]];
TD = [TD_with_time;[sim_time,TD_with_time(end,2)]];

t = 0:Ts:sim_time;
values = zeros(length(t),4);
values(1,:) = [FH_with_time(1,2),FC_with_time(1,2),FD_with_time(1,2),TD_with_time(1,2)];
% Merge all changing signal variable into one (ugly version).
for i=1:length(t)
% Main simulation.
    %Increase counters.
    actual_time = t(i);
    actual_FH = 0;
    actual_FC = 0;
    actual_FD = 0;
    actual_TD = 0;
    %Increase counters.
    if (actual_time >= FH(counter_FH,1)) && (counter_FH ~= FH_length)
        actual_FH = FH(counter_FH,2);
        counter_FH = counter_FH + 1;
    else
        actual_FH = FH(counter_FH-1,2);
    end

    if (actual_time >= FC(counter_FC,1)) && (counter_FC ~= FC_length)
        actual_FC = FC(counter_FC,2);
        counter_FC = counter_FC + 1;
    else
        actual_FC = FC(counter_FC-1,2);
    end

    if (actual_time >= FD(counter_FD,1)) && (counter_FD ~= FD_length)
        actual_FD = FD(counter_FD,2);
        counter_FD = counter_FD + 1;
    else
        actual_FD = FD(counter_FD-1,2);
    end

    if (actual_time >= TD(counter_TD,1)) && (counter_TD ~= TD_length)
        actual_TD = TD(counter_TD,2);
        counter_TD = counter_TD + 1;
    else
        actual_TD = TD(counter_TD-1,2);
    end

    values(i,:) = [actual_FH actual_FC actual_FD actual_TD];
    
    if actual_time == sim_time
        break;
    end
     %Find minimal time.
end
[y,t,x]=lsim(G, values, t ,x_0);
end

