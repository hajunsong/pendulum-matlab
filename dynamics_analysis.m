function dynamics_analysis

    global Y body h_temp g
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
%         body(i).err_pos = body(i).des_pos - body(i).qi;
%         body(i).err_pos_accum = body(i).err_pos_accum + body(i).err_pos*h_temp;
%         body(i).Tc_pos = body(i).Kp_pos*body(i).err_pos + body(i).Kd_pos*(body(i).err_pos - body(i).err_pos_prev)/h_temp + body(i).Ki_pos*body(i).err_pos_accum;
%         body(i).err_pos_prev = body(i).err_pos;
% 
%         body(i).err_vel = body(i).des_vel - body(i).dqi;
%         body(i).err_vel_accum = body(i).err_vel_accum + body(i).err_vel*h_temp;
%         body(i).Tc_vel = body(i).Kp_vel*body(i).err_vel + body(i).Kd_vel*(body(i).err_vel - body(i).err_vel_prev)/h_temp + body(i).Ki_vel*body(i).err_vel_accum;
%         body(i).err_vel_prev = body(i).err_vel;
% 
%         body(i).Tg = -body(i).mi*g*body(i).rhoip(2)*cos(body(i).qi);
% %         body(i).T_control = body(i).Tc_pos + body(i).Tc_vel + body(i).Tg;
%         body(i).T_control = 0;
    end
    
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
        Q(i,1) = body(i).Bi'*(body(i).Li - body(i).Ki*D_temp);
    end
    
    dYh = M\Q;
    for i = 1 : num_body
        body(i).ddqi = dYh(i);
    end
    
%     for i = 2 : num_body
%         %% System EQM
%         M = body(i).Bi'*body(i).Ki*body(i).Bi;
% %         if t_current >= 0.3 && t_current <= 0.32
% %             body(i).Td = -100;
% %         else
% %             body(i).Td = 0;
% %         end
%         body(i).Td = 0;
% 
%         Q = body(i).Bi'*(body(i).Li - body(i).Ki*body(i).Di) + body(i).T_control + body(i).Td;
% 
%         body(i).ddqi = M\Q;
%     
%         %% Collision Detect
%         body(i).r_hat_dot = body(i).T_control - body(i).Bi'*(body(i).Li - [body(i).Fic;body(i).rict*body(i).Fic] - body(i).Ki*body(i).Di) - body(i).r_hat;
%     end
    
    %% dqddq2Yp
    for i = 1 : num_body
        Yp(i,1) = body(i).dqi;
    end
    for i = 1 : num_body
        Yp(i+num_body,1) = body(i).ddqi;
    end

end