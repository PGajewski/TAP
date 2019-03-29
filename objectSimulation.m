function [t,x,y] = objectSimulation(FH_with_time,FC_with_time,FD_with_time, tau_H, tau_C, tau_D, tau_output, TH, TC, TD, sim_time, x_0)
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
assert((temp1(2) == 2 && temp2(2) == 2 && temp3(2) == 2),'Incompatible changing moment values!');
assert((isnumeric(tau_H) && isnumeric(tau_C) && isnumeric(tau_D) && isnumeric(tau_output)),'Transport delay must be a number!');
assert((isvector(x_0) && length(x_0) == 2),'Init condition must be a vector with size 2'); 

%Check last value change moment timestep.
assert((FH_with_time(end,1)<=sim_time) && (FC_with_time(end,1)<=sim_time) && (FC_with_time(end,1)<=sim_time),'One or more value changes is outside simulation range!');

%Add simulation end to changes moments.
FH = [FH_with_time;[sim_time,FH_with_time(end,2)]];
FC = [FC_with_time;[sim_time,FC_with_time(end,2)]];
FD = [FD_with_time;[sim_time,FD_with_time(end,2)]];

% Merge all changing signal variable into one (ugly version).
% Init.
FH_length = temp1(1)+1;
FC_length = temp2(1)+1;
FD_length = temp3(1)+1;
counter_C = 1;
counter_H = 1;
counter_D = 1;

values_with_time = [];
actual_time = 0;

while actual_time < sim_time
    %Find minimal time.
    actual_time = min([(FH(counter_H,1) + tau_H) (FC(counter_C,1) + tau_C) (FD(counter_D,1) + tau_D)]);
    
    %Add next timestamp with values.
    values_with_time = [values_with_time [actual_time;FH(counter_H,2);FC(counter_C,2);FD(counter_D,2)]];
        
    %Increase counters.
    if (actual_time == (FH(counter_H,1)+tau_H)) && (counter_H ~= FH_length)
        counter_H = counter_H + 1;
    end
    
    if (actual_time == (FC(counter_C,1)+tau_C)) && (counter_C ~= FC_length)
        counter_C = counter_C + 1;
    end
    
    if (actual_time == (FD(counter_D,1)+tau_D)) && (counter_D ~= FD_length)
        counter_D = counter_D + 1;
    end
end



% Main simulation.
x = [];
t = [];

% Simulate by all changes moments.
for i=1:(length(values_with_time)-1)
    stateHandler = @(t,x) stateFunction(t,x,values_with_time(2,i), values_with_time(3,i), values_with_time(4,i), TH, TC, TD);
    [temp1,temp2]=ode45(stateHandler,[values_with_time(1,i) values_with_time(1,i+1)],x_0);
    t = [t temp1(1:end-1,:)'];
    x = [x temp2(1:end-1,:)'];
    % Remember previous state for next step.
    x_0 = temp2(end,:);
end

process_length = length(t);

%Count output.
y = ones(1,process_length);
counter = 1;
for i=1:process_length
    if t(i) <= tau_output
       y(i) = x_0(2); 
    else
        %Count expected process time for get temperature.
        expected_time = t(i) - tau_output;

        %Found first value greater than expected.
       while counter < process_length
        	if t(counter) > expected_time
                y(i) = (expected_time - t(counter - 1)) * (x(2,counter) - x(2,counter - 1))/(t(counter) - t(counter - 1)) + x(2,counter - 1);
                break;
            elseif t(counter) == expected_time
                y(i) = x(2,counter);
                break;
            end
            counter = counter + 1;
        end
    end
end

end

