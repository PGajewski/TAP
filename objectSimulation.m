function [t,x,y] = objectSimulation(FH_with_time,FC_with_time,FD_with_time, TD_with_time, tau_H, tau_C, tau_D, tau_output, TH, TC, sim_time, x_0)
%objectSimulation Simulation of object
%   FH_with_time - matrix of changing FH values: new value (first column) and time step (second column).
%   FC_with_time - matrix of changing FC values: new value (first column) and time step (second column).
%   FD_with_time - matrix of changing FD values: new value (first column) and time step (second column).
%   TD_with_time - matrix of changing FD values: new value (first column) and time step (second column).
%   tau_H        - transport delay of signal FH.
%   tau_C        - transport delay of signal FC.
%   tau_D        - transport delay of signal FD.

% Check variables dimension.
temp1 = size(FH_with_time);
temp2 = size(FC_with_time);
temp3 = size(FD_with_time);
temp4 = size(TD_with_time);
assert((temp1(2) == 2 && temp2(2) == 2 && temp3(2) == 2 && temp4(2) == 2),'Incompatible changing moment values!');
assert((isnumeric(tau_H) && isnumeric(tau_C) && isnumeric(tau_D) && isnumeric(tau_output)),'Transport delay must be a number!');
assert((isvector(x_0) && length(x_0) == 2),'Init condition must be a vector with size 2'); 

%Check last value change moment timestep.
assert((FH_with_time(end,1)<=sim_time) && (FC_with_time(end,1)<=sim_time) && (FD_with_time(end,1)<=sim_time) && (TD_with_time(end,1)<=sim_time),'One or more value changes is outside simulation range!');

%Prepare ODE options.
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

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
counter_FC = 2;
counter_FH = 2;
counter_FD = 2;
counter_TD = 2;


values_with_time = [0;FH(1,2);FC(1,2);FD(1,2);TD(1,2)];
actual_time = min([FH(2,1) FC(2,1) FD(2,1) TD(2,1)]);

while 1
    %Add next timestamp with values.
    actual_FH = 0;
    actual_FC = 0;
    actual_FD = 0;
    actual_TD = 0;
    %Increase counters.
    if (actual_time ~= tau_H) && (actual_time >= (FH(counter_FH,1)+tau_H)) && (counter_FH ~= FH_length)
        actual_FH = FH(counter_FH,2);
        counter_FH = counter_FH + 1;
    else
        actual_FH = FH(counter_FH-1,2);
    end

    if (actual_time ~= tau_C) && (actual_time >= (FC(counter_FC,1)+tau_C)) && (counter_FC ~= FC_length)
        actual_FC = FC(counter_FC,2);
        counter_FC = counter_FC + 1;
    else
        actual_FC = FC(counter_FC-1,2);
    end

    if (actual_time ~= tau_D) && (actual_time >= (FD(counter_FD,1)+tau_D)) && (counter_FD ~= FD_length)
        actual_FD = FD(counter_FD,2);
        counter_FD = counter_FD + 1;
    else
        actual_FD = FD(counter_FD-1,2);
    end

    if (actual_time ~= tau_D) && (actual_time >= (TD(counter_TD,1)+tau_D)) && (counter_TD ~= TD_length)
        actual_TD = TD(counter_TD,2);
        counter_TD = counter_TD + 1;
    else
        actual_TD = TD(counter_TD-1,2);
    end

    values_with_time = [values_with_time [actual_time;actual_FH;actual_FC;actual_FD; actual_TD]];
    
    if actual_time == sim_time
        break;
    end
     %Find minimal time.
    actual_time = min([(FH(counter_FH,1) + tau_H) (FC(counter_FC,1) + tau_C) (FD(counter_FD,1) + tau_D) (TD(counter_TD,1) + tau_D)]);
end

% Main simulation.
x = [];
t = [];

s = size(values_with_time);
init_values = x_0;
% Simulate by all changes moments.
for i=1:(s(2)-1)
    stateHandler = @(t,x) stateFunction(t,x,values_with_time(2,i), values_with_time(3,i), values_with_time(4,i), TH, TC, values_with_time(5,i));
    [temp1,temp2]=ode45(stateHandler,[values_with_time(1,i) values_with_time(1,i+1)],init_values, options);
    t = [t temp1(1:end-1,:)'];
    x = [x temp2(1:end-1,:)'];
    % Remember previous state for next step.
    init_values = temp2(end,:);
end

t =[t sim_time];
x =[x init_values'];
process_length = length(t);

%Count output.
y = ones(2,process_length);
counter = 2;
%Temperature.
for i=1:process_length
    if t(i) < tau_output
        y(1,i)= x_0(1);
    else
        %Count expected process time for get temperature.
        expected_time = t(i) - tau_output;

        %Found first value greater than expected.
       while counter < process_length
            if t(counter) > expected_time
                y(1,i) = (expected_time - t(counter - 1)) * (x(1,counter) - x(1,counter - 1))/(t(counter) - t(counter - 1)) + x(1,counter - 1);
                break;
            elseif t(counter) == expected_time
                y(1,i) = x(1,counter);
                break;
            end
            counter = counter + 1;
       end
    end
end
%Hight.
y(2,:) = x(2,:);
end

