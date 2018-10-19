function dynamics_analysis

    global Y body g h
    global t_current num_body A0 C01 s01p
    global Yp
    
    %% Y2qdq
    for i = 1 : num_body
        body(i).qi = Y(i,1);
    end
    for i = 1 : num_body
         body(i).dqi = Y(i+num_body,1);
    end

    %% Body 1~n
    for i = 1 : num_body
        %% Orientation Body
        body(i).Aijpp = [cos(body(i).qi) -sin(body(i).qi) 0;
            sin(body(i).qi) cos(body(i).qi) 0;
            0 0 1];
        if i == 1
            A0 = [0,0,1;1,0,0;0,1,0];
            C01 = eye(3);
            body(i).Ai = A0*C01*body(i).Aijpp;
            body(i).Hi = A0*C01*[0;0;1];
        else
            body(i).Ai = body(i-1).Ai*body(i-1).Cij*body(i).Aijpp;
            body(i).Hi = body(i-1).Ai*body(i-1).Cij*[0;0;1];
        end

        %% Position Body
        if i == 1
            s01 = A0*s01p;
            body(i).ri = s01;
        else
            body(i-1).sij = body(i-1).Ai*body(i-1).sijp;
            body(i).ri = body(i-1).ri + body(i-1).sij;
        end
        body(i).rhoi = body(i).Ai*body(i).rhoip;
        body(i).ric = body(i).ri + body(i).rhoi;

        %% Velocity State Body
        body(i).rit = tilde(body(i).ri);
        body(i).Bi = [body(i).rit*body(i).Hi; body(i).Hi];
        if i == 1
            body(i).Yih = body(i).Bi*body(i).dqi;
        else
            body(i).Yih = body(i-1).Yih + body(i).Bi*body(i).dqi;
        end

        %% Cartesian Velocity Body
        body(i).Ti = [eye(3) -body(i).rit; zeros(3) eye(3)];
        body(i).Yib = body(i).Ti*body(i).Yih;
        body(i).dri = body(i).Yib(1:3,1);
        body(i).wi = body(i).Yib(4:6,1);
        body(i).wit = tilde(body(i).wi);
        body(i).dric = body(i).dri + body(i).wit*body(i).rhoi;

        %% Mass & Force Body
        body(i).Jic = body(i).Ai*body(i).Cii*body(i).Jip*(body(i).Ai*body(i).Cii)';
        body(i).rict = tilde(body(i).ric);
        body(i).drict = tilde(body(i).dric);

        body(i).Mih = [body(i).mi*eye(3) -body(i).mi*body(i).rict;
            body(i).mi*body(i).rict body(i).Jic-body(i).mi*body(i).rict*body(i).rict];

        body(i).Fic = [0;0;body(i).mi*g];
        body(i).Tic = [0;0;0];

        body(i).Qih = [body(i).Fic + body(i).mi*body(i).drict*body(i).wi;
            body(i).Tic + body(i).rict*body(i).Fic + body(i).mi*body(i).rict*body(i).drict*body(i).wi - body(i).wit*body(i).Jic*body(i).wi];

        %% Velocity Coupling Body
        body(i).drit = tilde(body(i).dri);
        if i == 1
            body(i).dHi = zeros(3,1);
        else
            body(i).dHi = body(i-1).wit*body(i).Hi;
        end
        body(i).Di = [body(i).drit*body(i).Hi + body(i).rit*body(i).dHi; body(i).dHi]*body(i).dqi;

        %% Control
        body(i).err_pos = body(i).des_pos - body(i).qi;
        body(i).err_pos_accum = body(i).err_pos_accum + body(i).err_pos*h;
        body(i).Tc_pos = body(i).Kp_pos*body(i).err_pos + body(i).Kd_pos*(body(i).err_pos - body(i).err_pos_prev)/h + body(i).Ki_pos*body(i).err_pos_accum;
        body(i).err_pos_prev = body(i).err_pos;

        body(i).err_vel = body(i).des_vel - body(i).dqi;
        body(i).err_vel_accum = body(i).err_vel_accum + body(i).err_vel*h;
        body(i).Tc_vel = body(i).Kp_vel*body(i).err_vel + body(i).Kd_vel*(body(i).err_vel - body(i).err_vel_prev)/h + body(i).Ki_vel*body(i).err_vel_accum;
        body(i).err_vel_prev = body(i).err_vel;

        body(i).T_control = body(i).Tc_pos + body(i).Tc_vel;
%         body(i).T_control = 0;
    end
%     
%     t1 = cos(body(1).qi);
%     t2 = sin(body(1).qi);
%     t3 = t1 * A0(3,1) + t2 * A0(3,2);
%     t4 = -t1 * A0(3,2) + t2 * A0(3,1);
%     t5 = t3 * body(1).rhoip(1) - t4 * body(1).rhoip(2) + A0(3,3) * body(1).rhoip(3) + body(1).ri(3);
%     t6 = t1 * A0(2,1) + t2 * A0(2,2);
%     t7 = -t1 * A0(2,2) + t2 * A0(2,1);
%     t8 = t6 * body(1).rhoip(1) - t7 * body(1).rhoip(2) + A0(2,3) * body(1).rhoip(3) + body(1).ri(2);
%     t9 = body(1).mi * t5;
%     t10 = body(1).mi * t8;
%     t11 = t1 * A0(1,1) + t2 * A0(1,2);
%     t1 = -t1 * A0(1,2) + t2 * A0(1,1);
%     t2 = -t1 * body(1).rhoip(2) + t11 * body(1).rhoip(1) + A0(1,3) * body(1).rhoip(3) + body(1).ri(1);
%     t12 = body(1).mi * t5;
%     t13 = body(1).mi * t2;
%     t14 = body(1).mi * t8;
%     t15 = body(1).mi * t2;
%     t16 = -t1 * body(1).Cii(2,1) + t11 * body(1).Cii(1,1) + A0(1,3) * body(1).Cii(3,1);
%     t17 = -t1 * body(1).Cii(2,2) + t11 * body(1).Cii(1,2) + A0(1,3) * body(1).Cii(3,2);
%     t1 = -t1 * body(1).Cii(2,3) + t11 * body(1).Cii(1,3) + A0(1,3) * body(1).Cii(3,3);
%     t11 = t1 * body(1).Jip(3,1) + t16 * body(1).Jip(1,1) + t17 * body(1).Jip(2,1);
%     t18 = t1 * body(1).Jip(3,2) + t16 * body(1).Jip(1,2) + t17 * body(1).Jip(2,2);
%     t19 = t1 * body(1).Jip(3,3) + t16 * body(1).Jip(1,3) + t17 * body(1).Jip(2,3);
%     t20 = t5 ^ 2;
%     t21 = t8 ^ 2;
%     t22 = t6 * body(1).Cii(1,1) - t7 * body(1).Cii(2,1) + A0(2,3) * body(1).Cii(3,1);
%     t23 = t6 * body(1).Cii(1,2) - t7 * body(1).Cii(2,2) + A0(2,3) * body(1).Cii(3,2);
%     t6 = t6 * body(1).Cii(1,3) - t7 * body(1).Cii(2,3) + A0(2,3) * body(1).Cii(3,3);
%     t7 = t3 * body(1).Cii(1,1) - t4 * body(1).Cii(2,1) + A0(3,3) * body(1).Cii(3,1);
%     t24 = t3 * body(1).Cii(1,2) - t4 * body(1).Cii(2,2) + A0(3,3) * body(1).Cii(3,2);
%     t3 = t3 * body(1).Cii(1,3) - t4 * body(1).Cii(2,3) + A0(3,3) * body(1).Cii(3,3);
%     t4 = t22 * body(1).Jip(1,1) + t23 * body(1).Jip(2,1) + t6 * body(1).Jip(3,1);
%     t25 = t22 * body(1).Jip(1,2) + t23 * body(1).Jip(2,2) + t6 * body(1).Jip(3,2);
%     t26 = t22 * body(1).Jip(1,3) + t23 * body(1).Jip(2,3) + t6 * body(1).Jip(3,3);
%     t27 = t2 ^ 2;
%     t28 = t24 * body(1).Jip(2,1) + t3 * body(1).Jip(3,1) + t7 * body(1).Jip(1,1);
%     t29 = t24 * body(1).Jip(2,2) + t3 * body(1).Jip(3,2) + t7 * body(1).Jip(1,2);
%     t30 = t24 * body(1).Jip(2,3) + t3 * body(1).Jip(3,3) + t7 * body(1).Jip(1,3);
%     body(1).Mih = [body(1).mi 0 0 0 t9 -t10; 
%         0 body(1).mi 0 -t12 0 t13; 
%         0 0 body(1).mi t14 -t15 0; 
%         0 -t9 t10 t1 * t19 + t11 * t16 + t17 * t18 + (t20 + t21) * body(1).mi -t10 * t2 + t11 * t22 + t18 * t23 + t19 * t6 t11 * t7 + t18 * t24 + t19 * t3 - t9 * t2; 
%         t12 0 -t13 t1 * t26 - t13 * t8 + t16 * t4 + t17 * t25 t22 * t4 + t23 * t25 + t26 * t6 + (t20 + t27) * body(1).mi -t10 * t5 + t24 * t25 + t26 * t3 + t4 * t7; 
%         -t14 t15 0 t1 * t30 - t15 * t5 + t16 * t28 + t17 * t29 t22 * t28 + t23 * t29 + t30 * t6 - t9 * t8 t24 * t29 + t28 * t7 + t3 * t30 + (t21 + t27) * body(1).mi;];
% 
%     t1 = cos(body(1).qi);
%     t2 = sin(body(1).qi);
%     t3 = t1 * A0(1,1) + t2 * A0(1,2);
%     t4 = -t1 * A0(1,2) + t2 * A0(1,1);
%     t5 = t3 * body(1).rhoip(1);
%     t6 = t4 * body(1).rhoip(2);
%     t7 = A0(1,3) * body(1).rhoip(3);
%     t8 = t1 * A0(2,1) + t2 * A0(2,2);
%     t9 = -t1 * A0(2,2) + t2 * A0(2,1);
%     t10 = t8 * body(1).rhoip(1);
%     t11 = t9 * body(1).rhoip(2);
%     t12 = A0(2,3) * body(1).rhoip(3);
%     t13 = t5 - t6 + body(1).ri(1) + t7;
%     t14 = -t10 + t11 - t12 - body(1).ri(2);
%     t15 = body(i).dqi * (t13 * A0(2,3) + t14 * A0(1,3) + A0(1,3) * body(1).ri(2) - A0(2,3) * body(1).ri(1));
%     t16 = t1 * A0(3,1) + t2 * A0(3,2);
%     t1 = -t1 * A0(3,2) + t2 * A0(3,1);
%     t2 = t16 * body(1).rhoip(1);
%     t17 = t1 * body(1).rhoip(2);
%     t18 = A0(3,3) * body(1).rhoip(3);
%     t19 = t2 - t17 + t18;
%     t13 = body(i).dqi * (t13 * A0(3,3) + (-t19 - body(1).ri(3)) * A0(1,3) + A0(1,3) * body(1).ri(3) - A0(3,3) * body(1).ri(1));
%     t20 = t13 * A0(3,3);
%     t21 = body(1).mi * body(i).dqi;
%     t14 = body(i).dqi * (t14 * A0(3,3) + (t19 + body(1).ri(3)) * A0(2,3) - A0(2,3) * body(1).ri(3) + A0(3,3) * body(1).ri(2));
%     t19 = t14 * A0(3,3);
%     t22 = A0(2,3) * t14;
%     t23 = A0(1,3) * t13;
%     t10 = t10 + body(1).ri(2) - t11 + t12;
%     t2 = body(1).ri(3) + t2 - t17 + t18;
%     t11 = t8 * body(1).Cii(1,1) - t9 * body(1).Cii(2,1) + A0(2,3) * body(1).Cii(3,1);
%     t12 = t8 * body(1).Cii(1,2) - t9 * body(1).Cii(2,2) + A0(2,3) * body(1).Cii(3,2);
%     t8 = t8 * body(1).Cii(1,3) - t9 * body(1).Cii(2,3) + A0(2,3) * body(1).Cii(3,3);
%     t9 = t11 * body(1).Jip(1,1) + t12 * body(1).Jip(2,1) + t8 * body(1).Jip(3,1);
%     t17 = t3 * body(1).Cii(1,1) - t4 * body(1).Cii(2,1) + A0(1,3) * body(1).Cii(3,1);
%     t18 = t11 * body(1).Jip(1,2) + t12 * body(1).Jip(2,2) + t8 * body(1).Jip(3,2);
%     t24 = t3 * body(1).Cii(1,2) - t4 * body(1).Cii(2,2) + A0(1,3) * body(1).Cii(3,2);
%     t25 = t11 * body(1).Jip(1,3) + t12 * body(1).Jip(2,3) + t8 * body(1).Jip(3,3);
%     t3 = t3 * body(1).Cii(1,3) - t4 * body(1).Cii(2,3) + A0(1,3) * body(1).Cii(3,3);
%     t4 = t17 * t9 + t18 * t24 + t25 * t3;
%     t26 = -t1 * body(1).Cii(2,1) + t16 * body(1).Cii(1,1) + A0(3,3) * body(1).Cii(3,1);
%     t27 = -t1 * body(1).Cii(2,2) + t16 * body(1).Cii(1,2) + A0(3,3) * body(1).Cii(3,2);
%     t1 = -t1 * body(1).Cii(2,3) + t16 * body(1).Cii(1,3) + A0(3,3) * body(1).Cii(3,3);
%     t16 = t1 * body(1).Jip(3,1) + t26 * body(1).Jip(1,1) + t27 * body(1).Jip(2,1);
%     t28 = t1 * body(1).Jip(3,2) + t26 * body(1).Jip(1,2) + t27 * body(1).Jip(2,2);
%     t29 = t1 * body(1).Jip(3,3) + t26 * body(1).Jip(1,3) + t27 * body(1).Jip(2,3);
%     t30 = t16 * t17 + t24 * t28 + t29 * t3;
%     t31 = t11 * t9 + t12 * t18 + t25 * t8;
%     t32 = t11 * t16 + t12 * t28 + t29 * t8;
%     t9 = t1 * t25 + t18 * t27 + t26 * t9;
%     t16 = t1 * t29 + t16 * t26 + t27 * t28;
%     t5 = t5 - t6 + body(1).ri(1) + t7;
%     t6 = t17 * body(1).Jip(1,1) + t24 * body(1).Jip(2,1) + t3 * body(1).Jip(3,1);
%     t7 = t17 * body(1).Jip(1,2) + t24 * body(1).Jip(2,2) + t3 * body(1).Jip(3,2);
%     t18 = t17 * body(1).Jip(1,3) + t24 * body(1).Jip(2,3) + t3 * body(1).Jip(3,3);
%     t3 = t17 * t6 + t18 * t3 + t24 * t7;
%     t8 = t11 * t6 + t12 * t7 + t18 * t8;
%     t1 = t1 * t18 + t26 * t6 + t27 * t7;
%     t6 = [-body(i).dqi * ((body(1).mi * t15 * t5 - body(i).dqi * (t3 * A0(2,3) - t4 * A0(1,3))) * A0(1,3) + (body(1).mi * t15 * t10 - body(i).dqi * (-t31 * A0(1,3) + t8 * A0(2,3))) * A0(2,3) + (body(1).mi * (t10 * t13 + t14 * t5) - body(i).dqi * (t1 * A0(2,3) - t9 * A0(1,3))) * A0(3,3))];
%     body(i).Qih = [[t21 * (t15 * A0(2,3) + t20)] [-t21 * (t15 * A0(1,3) + t19)] [body(1).mi * ((-t23 + t22) * body(i).dqi + g)] [body(i).dqi * (body(i).dqi * (t31 * A0(3,3) - t32 * A0(2,3)) * A0(2,3) + body(i).dqi * (-t16 * A0(2,3) + t9 * A0(3,3)) * A0(3,3)) + (t10 * (t22 * body(i).dqi + g) + t19 * body(i).dqi * t2) * body(1).mi + body(i).dqi * (body(1).mi * (-t10 * t13 + t15 * t2) + body(i).dqi * (-t30 * A0(2,3) + t4 * A0(3,3))) * A0(1,3)] [-body(i).dqi * (body(i).dqi * (t3 * A0(3,3) - t30 * A0(1,3)) * A0(1,3) + body(i).dqi * (t1 * A0(3,3) - t16 * A0(1,3)) * A0(3,3)) - (t5 * (-t23 * body(i).dqi + g) - t20 * body(i).dqi * t2) * body(1).mi - body(i).dqi * (-body(1).mi * (-t14 * t5 + t15 * t2) + body(i).dqi * (-t32 * A0(1,3) + t8 * A0(3,3))) * A0(2,3)] t6;]';

    
    %% system EQM
    for i = num_body : -1 : 1
        body(i).Ki = body(i).Mih;
        body(i).Li = body(i).Qih;
        if i ~= num_body
            body(i).Ki = body(i).Ki + body(i+1).Ki;
            body(i).Li = body(i).Li + body(i+1).Li - body(i+1).Ki*body(i+1).Di;
        end
    end
    
    for i = 1 : num_body
        for j = 1 : num_body
            if i == j
                M(i, j) = body(i).Bi'*body(i).Ki*body(i).Bi;
            elseif i < j
                M(i, j) = body(i).Bi'*body(j).Ki*body(j).Bi;
            elseif i > j
                M(i, j) = body(i).Bi'*body(i).Ki*body(j).Bi;
            end
        end
        
        D_temp = zeros(6,1);
        for j = 1 : i
            D_temp = D_temp + body(j).Di;
        end
        % gravity compensation force
        body(i).Fg = -[body(i).Fic;body(i).rict*body(i).Fic];
        body(i).Tg = body(i).Bi'*(body(i).Fg - body(i).Ki*D_temp);
        
        body(i).T_in = body(i).T_control;
        
        body(i).Ta = body(i).Tg + body(i).T_in;
%         body(i).Ta = 0;
        
        if t_current >= 2 && t_current <= 2.05
            body(i).Td = -15;
        else
            body(i).Td = 0;
        end
        body(i).Td = 0;
        
        Q(i,1) = body(i).Bi'*(body(i).Li - body(i).Ki*D_temp) + body(i).Ta + body(i).Td;
        
        body(i).yp = body(i).Ta - body(i).Tg - body(i).r_hat ;
    end
%     
%     t1 = (A0(3,1) ^ 2);
%     t2 = (A0(3,2) ^ 2);
%     t3 = (A0(2,1) ^ 2);
%     t4 = (A0(2,2) ^ 2);
%     t5 = (A0(1,2) * A0(3,2));
%     t6 = (A0(1,1) * A0(3,1));
%     t7 = (A0(1,2) ^ 2);
%     t8 = (A0(1,1) ^ 2);
%     t9 = (A0(1,3) ^ 2);
%     t10 = (A0(3,3) ^ 2);
%     t11 = (A0(2,3) ^ 2);
%     t12 = 2 * A0(3,3);
%     t13 = -t12 * (t6 - t5) * A0(1,3) - (-t4 + t3 - t1 + t2) * t9 - (-t8 + t7 - t4 + t3) * t10 - t11 * (-t8 + t7 - t1 + t2);
%     t14 = A0(2,1) * A0(2,2);
%     t15 = A0(3,1) * A0(3,2);
%     t16 = A0(1,1) * A0(1,2);
%     t14 = -t11 * (t16 + t15) + t9 * (t14 - t15) + (A0(1,3) * (A0(1,1) * A0(3,2) + A0(1,2) * A0(3,1)) + (-t16 + t14) * A0(3,3)) * A0(3,3);
%     t15 = body(1).rhoip(1) ^ 2;
%     t16 = body(1).rhoip(2) ^ 2;
%     t17 = t14 * body(1).rhoip(2);
%     t18 = cos(body(1).qi);
%     t19 = sin(body(1).qi);
%     t20 = A0(2,3) * body(1).rhoip(3);
%     t21 = body(1).ri(2) + t20;
%     t22 = body(1).ri(2) - t20;
%     t23 = t21 * (t9 * A0(2,1) + A0(3,3) * (A0(2,1) * A0(3,3) - A0(2,3) * A0(3,1))) + A0(1,1) * A0(2,3) * t22 * A0(1,3);
%     t21 = t21 * (t9 * A0(2,2) + A0(3,3) * (A0(2,2) * A0(3,3) - A0(2,3) * A0(3,2))) + A0(1,2) * A0(2,3) * t22 * A0(1,3);
%     M = ((-4 * t17 * body(1).rhoip(1) - t13 * (-t15 + t16)) * t18 ^ 2 - (4 * t20 * body(1).ri(2) * t9) + (2 * t17 * body(1).rhoip(1)) + (t15 * (-t5 * t12 * A0(1,3) - t10 * (-t7 + t4) + t11 * (t7 + t2) - t9 * (t4 - t2))) + (t16 * (-t6 * t12 * A0(1,3) - t10 * (-t8 + t3) + t11 * (t8 + t1) - t9 * (t3 - t1)))) * body(1).mi + (t11 * body(1).Jic(2,2)) + ((body(1).Jic(1,1) * A0(1,3) + (body(1).Jic(1,3) + body(1).Jic(3,1)) * A0(3,3) + A0(2,3) * (body(1).Jic(1,2) + body(1).Jic(2,1))) * A0(1,3)) + ((A0(2,3) * (body(1).Jic(2,3) + body(1).Jic(3,2)) + body(1).Jic(3,3) * A0(3,3)) * A0(3,3)) - 0.2e1 * body(1).mi * (t18 * ((t23 * body(1).rhoip(1)) + (t21 * body(1).rhoip(2)) + (t14 * t15 + (t13 * body(1).rhoip(1) - t17) * body(1).rhoip(2)) * t19) + t19 * (t21 * body(1).rhoip(1) - t23 * body(1).rhoip(2)));
% 
%     t1 = (A0(2,1) * A0(3,1) - A0(2,2) * A0(3,2));
%     t2 = (A0(3,2) * A0(2,1) + A0(2,2) * A0(3,1));
%     t3 = (body(1).rhoip(1) ^ 2);
%     t4 = (body(1).rhoip(2) ^ 2);
%     t5 = (t3 - t4);
%     t6 = 2;
%     t7 = t6 * body(1).rhoip(2);
%     t8 = -t7 * t2 * body(1).rhoip(1) - t1 * t5;
%     t9 = (A0(1,1) * A0(2,1) - A0(1,2) * A0(2,2));
%     t10 = (A0(2,2) * A0(1,1) + A0(1,2) * A0(2,1));
%     t11 = -t7 * t10 * body(1).rhoip(1) - t9 * t5;
%     t12 = (A0(2,2) + A0(2,1));
%     t13 = (A0(2,1) - A0(2,2));
%     t14 = (A0(1,1) * A0(3,1) - A0(1,2) * A0(3,2));
%     t15 = (A0(1,1) * A0(3,2) + A0(1,2) * A0(3,1));
%     t16 = (A0(1,3) ^ 2);
%     t17 = (A0(1,3) * t16);
%     t18 = (A0(3,3) ^ 2);
%     t19 = (A0(2,3) ^ 2);
%     t20 = (A0(2,3) * t19);
%     t21 = cos(body(1).qi);
%     t1 = ((t2 * t5) / 0.2e1 - (body(1).rhoip(2) * t1 * body(1).rhoip(1)));
%     t2 = -t7 * t9 * body(1).rhoip(1) + t10 * t5;
%     t9 = A0(2,2) * body(1).rhoip(1);
%     t10 = A0(2,1) * body(1).rhoip(2);
%     t22 = t10 - t9;
%     t23 = A0(2,1) * body(1).rhoip(1) + A0(2,2) * body(1).rhoip(2);
%     t24 = A0(2,3) * t22;
%     t25 = sin(body(1).qi);
%     t26 = A0(3,1) * body(1).rhoip(1) + A0(3,2) * body(1).rhoip(2);
%     t27 = A0(2,3) * body(1).rhoip(3);
%     t28 = body(1).ri(2) + t27;
%     t29 = A0(1,1) * body(1).rhoip(1) + A0(1,2) * body(1).rhoip(2);
%     t19 = t18 * t28 - t19 * (body(1).ri(2) - t27);
%     t30 = (A0(2,3) - A0(3,3)) * (A0(2,3) + A0(3,3)) - t16;
%     t31 = body(1).dqi ^ 2;
%     t32 = A0(3,1) * body(1).rhoip(2) - A0(3,2) * body(1).rhoip(1);
%     t33 = A0(1,1) * body(1).rhoip(2) - A0(1,2) * body(1).rhoip(1);
%     t34 = t22 * t32;
%     t1 = (body(1).mi * (t21 * (-0.4e1 * t31 * ((t1 * t17) + (-(t18 * t2) / 0.2e1 + (A0(1,3) * (t1 * A0(3,3) - t24 * t23))) * A0(3,3) + (t20 * (-t7 * t14 * body(1).rhoip(1) + t15 * t5)) / 0.2e1 - (A0(3,3) * t2 * t16) / 0.2e1) * t25 + ((t23 * A0(1,3) - t29 * A0(2,3)) * g) + (-t6 * t28 * (t30 * t29 * A0(3,3) + t17 * t26) + (-4 * A0(2,3) * t28 * t23 * A0(3,3) - t6 * t26 * t19) * A0(1,3)) * t31) + t25 * (-(g * (t22 * A0(1,3) - t33 * A0(2,3))) + t31 * (t6 * t28 * (t30 * t33 * A0(3,3) + t17 * t32) + (t6 * t32 * t19 + 4 * t24 * t28 * A0(3,3)) * A0(1,3)))));
%     Q = t6 * ((t31 * (t8 * t17 + (-t11 * t18 + A0(1,3) * (-t11 * A0(1,3) + t8 * A0(3,3) - (t12 * body(1).rhoip(1) - t13 * body(1).rhoip(2)) * (t12 * body(1).rhoip(2) + t13 * body(1).rhoip(1)) * A0(2,3))) * A0(3,3) - (t7 * t15 * body(1).rhoip(1) + t14 * t5) * t20) * t21 ^ 2 + t31 * (t33 * (t18 * t22 * A0(3,3) - t20 * t32) + ((t22 * A0(3,3) * t33 - t34 * A0(1,3)) * A0(1,3) - A0(3,3) * ((t3 * A0(2,2) ^ 2 + t4 * A0(2,1) ^ 2 + t6 * (-t10 * t9 + t27 * t28)) * A0(2,3) + t34 * A0(3,3))) * A0(1,3)) + (g * body(1).ri(2) * A0(1,3))) * body(1).mi - A0(2,3) * A0(3,3) * t31 * (body(1).Jic(1,1) * A0(1,3) + A0(2,3) * body(1).Jic(1,2) + A0(3,3) * body(1).Jic(1,3))) + t1;
%     
    dYh = M\Q;
    for i = 1 : num_body
        body(i).ddqi = dYh(i);
    end
    
    %% dqddq2Yp
    for i = 1 : num_body
        Yp(i,1) = body(i).dqi;
    end
    for i = 1 : num_body
        Yp(i+num_body,1) = body(i).ddqi;
    end
    for i = 1 : num_body
        Yp(num_body*2 + i, 1) = body(i).yp;
    end

end