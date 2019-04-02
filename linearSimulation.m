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

%Add simulation end to changes moments.
FH = [FH_with_time;[sim_time,FH_with_time(end,2)]];
FC = [FC_with_time;[sim_time,FC_with_time(end,2)]];
FD = [FD_with_time;[sim_time,FD_with_time(end,2)]];
TD = [TD_with_time;[sim_time,TD_with_time(end,2)]];

% Merge all changing signal variable into one (ugly version).
% Init.
FH_length = temp1(1)+1;
FC_length = temp2(1)+1;
FD_length = temp3(1)+1;
TD_length = temp4(1)+1;
counter_FC = 1;
counter_FH = 1;
counter_FD = 1;
counter_TD = 1;


values_with_time = [];

for actual_time=0:Ts:sim_time
    values_with_time = [values_with_time [actual_time;FH(counter_FH,2);FC(counter_FC,2);FD(counter_FD,2); TD(counter_FD,2)]];

    %Increase counters.
    if (actual_time == FH(counter_FH,1)) && (counter_FH ~= FH_length)
        counter_FH = counter_FH + 1;
    end
    
    if (actual_time == FC(counter_FC,1)) && (counter_FC ~= FC_length)
        counter_FC = counter_FC + 1;
    end
    
    if (actual_time == FD(counter_FD,1)) && (counter_FD ~= FD_length)
        counter_FD = counter_FD + 1;
    end
    
    if (actual_time == TD(counter_TD,1)) && (counter_TD ~= TD_length)
        counter_TD = counter_TD + 1;
    end
end

% Simulation.
[y,t,x] = lsim(G,values_with_time(2:5,:)',values_with_time(1,:)','zoh');

end

