function [t,values] = countErrorFunction(t1,y1,t2,y2)
%countErrorFunction Function to plot diference between two vectors.
%   Detailed explanation goes here
temp1 = size(t1);
temp2 = size(t2);
temp3 = size(y1);
temp4 = size(y2);

assert((temp1(2) == temp3(2)) && (temp2(2) == temp4(2)),'Incompatible vectors length!');
t = [];
values = [];
counter_y1 = 1;
counter_y2 = 1;
while 1
    if t1(1,counter_y1) == t2(1,counter_y2)
        t = [t; t1(1,counter_y1)];
        values = [values; y2(1,counter_y2) - y1(1,counter_y1)];
        if counter_y1 == temp2(2) || counter_y2 == temp4(2)
            break;
        end
        counter_y1 = counter_y1 + 1;
        counter_y2 = counter_y2 + 1;
    elseif t1(counter_y1) > t2(counter_y2)
        if counter_y2 == temp4(2)
            break;
        end
        counter_y2 = counter_y2 + 1;  
    else
       if counter_y1 == temp2(2)
            break;
       end
        counter_y1 = counter_y1 + 1; 
    end
end
end

