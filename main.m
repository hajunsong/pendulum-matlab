clc; clear all; close all;

global t_current h_temp start_time end_time h Y Yp body body1 body2 body3 body4 body5 body6 num_body

format long g

read_data;

num_body = 1;
switch(num_body)
    case 1
        body = body1;
    case 2
        body = [body1, body2];
    case 3
        body = [body1, body2, body3];
    case 4
        body = [body1, body2, body3, body4];
    case 5
        body = [body1, body2, body3, body4, body5];
    case 6
        body = [body1, body2, body3, body4, body5, body6];
end

define_Y_vector;

intcount = 1;
t_current = 0;
indx = 1;
h_temp = h/4;

delta_theta = body(1).des_pos - body1.qi;
a = 6*delta_theta/3^5;
b = -15*delta_theta/3^4;
c = 10*delta_theta/3^3;

t_inter = 2;
flag = 0;

file_name = sprintf('matlab_body%d.txt', num_body);
fp = fopen(file_name, 'w+');

data = zeros(1, num_body);
while(t_current <= end_time)
%     if t_current <= t_inter
%         body(1).des_pos = a*t_current^5 + b*t_current^4 + c*t_current^3;
%     end
    
    dynamics_analysis;
    
    [Y_next, t_next, intcount] = absh3(t_current, Y, Yp, h, intcount);
    
    data(indx,1) = t_current;
    fprintf(fp, '%.5f\t',t_current);
    for i = 1 : num_body
        data(indx, 3*i-1) = body(i).qi;
        data(indx, 3*i) = body(i).dqi;
        data(indx, 3*i+1) = body(i).ddqi;
        data(indx, 3*i+2) = body(i).r_hat;
        fprintf(fp, '%.5f\t%.5f\t%.5f\t', body(i).qi, body(i).dqi, body(i).ddqi);
    end
    fprintf(fp, '\n');
   
    u_vec = [1;0;0];
    body(1).r_hat = body(1).K*(Y_next(2*num_body + 1,1) - u_vec'*body(1).Jic*u_vec*0.5*body(1).dqi*body(1).dqi);
    
%     if t_current > t_inter && flag == 0
%         body(1).des_pos = body(1).qi;
%         flag = 1;
%     end
    
    Y = Y_next;
    h_temp = t_next - t_current;
    t_current = t_next
    indx = indx + 1;
end
fclose(fp);
save('matlab_save_data','data')

plotting;
% plotting2;