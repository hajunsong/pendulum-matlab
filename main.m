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

delta_w = body(1).des_vel;
t1 = 0.5;
a = -2*delta_w/t1^3;
b = 3*delta_w/t1^2;

t_inter = 2;
flag = 0;

Y_old = zeros(size(Y,1),1);
Yp_old = zeros(size(Y,1),1);

file_name = sprintf('matlab_body%d.txt', num_body);
fp = fopen(file_name, 'w+');

data = zeros(1, num_body);
while(t_current <= end_time)
    if t_current <= t1
        body(1).des_vel = a*t_current^3 + b*t_current^2;
    end
    
    dynamics_analysis;
    
%     [Y, t_current, intcount] = absh3(t_current, Y, Yp, h, intcount);
    Y = Y_old + Yp_old*h + 0.5*h*(Yp - Yp_old);
    
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
   
    body(1).p = 0.5*body(1).dqi^2*(body(1).Jic(1,1)*body(1).wi(1)^2 + body(1).Jic(2,2)*body(1).wi(2)^2 + body(1).Jic(3,3)*body(1).wi(3)^2);
    body(1).r_hat = body(1).K*(Y(2*num_body + 1,1) - body(1).p);
    
%     if t_current > t_inter && flag == 0
%         body(1).des_pos = body(1).qi;
%         flag = 1;
%     end
    
%     Y = Y_next;
%     h_temp = t_next - t_current;
%     t_current = t_next
    Y_old = Y;
    Yp_old = Yp;
    t_current = t_current + h
    indx = indx + 1;
end
fclose(fp);
save('matlab_save_data','data')

plotting;
% plotting2;