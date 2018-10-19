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
    end
    
    
    %% Orientation Body
    body(1).Aijpp = [cos(body(1).qi) -sin(body(1).qi) 0; sin(body(1).qi) cos(body(1).qi) 0; 0 0 1];
    body(2).Aijpp = [cos(body(2).qi) -sin(body(2).qi) 0; sin(body(2).qi) cos(body(2).qi) 0; 0 0 1];
    body(3).Aijpp = [cos(body(3).qi) -sin(body(3).qi) 0; sin(body(3).qi) cos(body(3).qi) 0; 0 0 1];
    
    A0 = [0,0,1;1,0,0;0,1,0];
    C01 = eye(3);
    u_vec = [0;0;1];
    
    A0_C01 = A0*C01;
    body(1).Ai = A0_C01*body(1).Aijpp;
    body(1).Hi = A0_C01*u_vec;
    
    A1_C12 = body(1).Ai*body(1).Cij;
    body(2).Ai = A1_C12*body(2).Aijpp;
    body(2).Hi = A1_C12*u_vec;
    
    A2_C23 = body(2).Ai*body(2).Cij;
    body(3).Ai = A2_C23*body(3).Aijpp;
    body(3).Hi = A2_C23*u_vec;
    
    %% Position Body
    s01 = A0*s01p;
    body(1).sij = body(1).Ai*body(1).sijp;
    body(2).sij = body(2).Ai*body(2).sijp;
    
    body(1).ri = s01;
    body(2).ri = body(1).ri + body(1).sij;
    body(3).ri = body(2).ri + body(2).sij;
    
    body(1).rhoi = body(1).Ai*body(1).rhoip;
    body(1).ric = body(1).ri + body(1).rhoi;
    body(2).rhoi = body(2).Ai*body(2).rhoip;
    body(2).ric = body(2).ri + body(2).rhoi;
    body(3).rhoi = body(3).Ai*body(3).rhoip;
    body(3).ric = body(3).ri + body(3).rhoi;
    
    %% Velocity State Body
    body(1).rit = tilde(body(1).ri);
    body(2).rit = tilde(body(2).ri);
    body(3).rit = tilde(body(3).ri);
    
    body(1).Bi = [body(1).rit*body(1).Hi; body(1).Hi];
    body(2).Bi = [body(2).rit*body(2).Hi; body(2).Hi];
    body(3).Bi = [body(3).rit*body(3).Hi; body(3).Hi];
    
    body(1).Yih =               body(1).Bi*body(1).dqi;
    body(2).Yih = body(1).Yih + body(2).Bi*body(2).dqi;
    body(3).Yih = body(2).Yih + body(3).Bi*body(3).dqi;
    
    %% Cartesian Velocity Body
    body(1).Ti = [eye(3) -body(1).rit; zeros(3) eye(3)];
    body(1).Yib = body(1).Ti*body(1).Yih;
    body(1).dri = body(1).Yib(1:3,1);
    body(1).wi = body(1).Yib(4:6,1);
    body(1).wit = tilde(body(1).wi);
    body(1).dric = body(1).dri + body(1).wit*body(1).rhoi;
    
    body(2).Ti = [eye(3) -body(2).rit; zeros(3) eye(3)];
    body(2).Yib = body(2).Ti*body(2).Yih;
    body(2).dri = body(2).Yib(1:3,1);
    body(2).wi = body(2).Yib(4:6,1);
    body(2).wit = tilde(body(2).wi);
    body(2).dric = body(2).dri + body(2).wit*body(2).rhoi;
    
    body(3).Ti = [eye(3) -body(3).rit; zeros(3) eye(3)];
    body(3).Yib = body(3).Ti*body(3).Yih;
    body(3).dri = body(3).Yib(1:3,1);
    body(3).wi = body(3).Yib(4:6,1);
    body(3).wit = tilde(body(3).wi);
    body(3).dric = body(3).dri + body(3).wit*body(3).rhoi;
    
    %% Mass & Force Body
    body(1).Jic = body(1).Ai*body(1).Cii*body(1).Jip*(body(1).Ai*body(1).Cii)';
    body(1).rict = tilde(body(1).ric);
    
    body(1).Mih = [body(1).mi*eye(3) -body(1).mi*body(1).rict;
        body(1).mi*body(1).rict body(1).Jic-body(1).mi*body(1).rict*body(1).rict];
    body(1).Fic = [0;0;body(1).mi*g];
    body(1).Fg = -[body(1).Fic;body(1).rict*body(1).Fic];
    
    body(2).Jic = body(2).Ai*body(2).Cii*body(2).Jip*(body(2).Ai*body(2).Cii)';
    body(2).rict = tilde(body(2).ric);
    
    body(2).Mih = [body(2).mi*eye(3) -body(2).mi*body(2).rict;
        body(2).mi*body(2).rict body(2).Jic-body(2).mi*body(2).rict*body(2).rict];
    body(2).Fic = [0;0;body(2).mi*g];
    body(2).Fg = -[body(2).Fic;body(2).rict*body(2).Fic];
    
    body(3).Jic = body(3).Ai*body(3).Cii*body(3).Jip*(body(3).Ai*body(3).Cii)';
    body(3).rict = tilde(body(3).ric);
    
    body(3).Mih = [body(3).mi*eye(3) -body(3).mi*body(3).rict;
        body(3).mi*body(3).rict body(3).Jic-body(3).mi*body(3).rict*body(3).rict];
    body(3).Fic = [0;0;body(3).mi*g];
    body(3).Fg = -[body(3).Fic;body(3).rict*body(3).Fic];
	
	%% Velocity Coupling Body
	body(1).drit = tilde(body(1).dri);
	body(2).drit = tilde(body(2).dri);
	body(3).drit = tilde(body(3).dri);
	
	body(1).dHi = zeros(3,1);
	body(2).dHi = body(1).wit*body(2).Hi;
	body(3).dHi = body(2).wit*body(3).Hi;

	body(1).Di = [body(1).drit*body(1).Hi + body(1).rit*body(1).dHi; body(1).dHi]*body(1).dqi;
	body(2).Di = [body(2).drit*body(2).Hi + body(2).rit*body(2).dHi; body(2).dHi]*body(2).dqi;
	body(3).Di = [body(3).drit*body(3).Hi + body(3).rit*body(3).dHi; body(3).dHi]*body(3).dqi;
	
	%% Control
	body(1).err_pos = body(1).des_pos - body(1).qi;
	body(1).err_pos_accum = body(1).err_pos_accum + body(1).err_pos*h;
	body(1).Tc_pos = body(1).Kp_pos*body(1).err_pos + body(1).Kd_pos*(body(1).err_pos - body(1).err_pos_prev)/h + body(1).Ki_pos*body(1).err_pos_accum;
	body(1).err_pos_prev = body(1).err_pos;

	body(1).err_vel = body(1).des_vel - body(1).dqi;
	body(1).err_vel_accum = body(1).err_vel_accum + body(1).err_vel*h;
	body(1).Tc_vel = body(1).Kp_vel*body(1).err_vel + body(1).Kd_vel*(body(1).err_vel - body(1).err_vel_prev)/h + body(1).Ki_vel*body(1).err_vel_accum;
	body(1).err_vel_prev = body(1).err_vel;

	body(1).T_control = body(1).Tc_pos + body(1).Tc_vel;
    body(1).T_control = 0;

	body(2).err_pos = body(2).des_pos - body(2).qi;
	body(2).err_pos_accum = body(2).err_pos_accum + body(2).err_pos*h;
	body(2).Tc_pos = body(2).Kp_pos*body(2).err_pos + body(2).Kd_pos*(body(2).err_pos - body(2).err_pos_prev)/h + body(2).Ki_pos*body(2).err_pos_accum;
	body(2).err_pos_prev = body(2).err_pos;

	body(2).err_vel = body(2).des_vel - body(2).dqi;
	body(2).err_vel_accum = body(2).err_vel_accum + body(2).err_vel*h;
	body(2).Tc_vel = body(2).Kp_vel*body(2).err_vel + body(2).Kd_vel*(body(2).err_vel - body(2).err_vel_prev)/h + body(2).Ki_vel*body(2).err_vel_accum;
	body(2).err_vel_prev = body(2).err_vel;

	body(2).T_control = body(2).Tc_pos + body(2).Tc_vel;
    body(2).T_control = 0;

	body(3).err_pos = body(3).des_pos - body(3).qi;
	body(3).err_pos_accum = body(3).err_pos_accum + body(3).err_pos*h;
	body(3).Tc_pos = body(3).Kp_pos*body(3).err_pos + body(3).Kd_pos*(body(3).err_pos - body(3).err_pos_prev)/h + body(3).Ki_pos*body(3).err_pos_accum;
	body(3).err_pos_prev = body(3).err_pos;

	body(3).err_vel = body(3).des_vel - body(3).dqi;
	body(3).err_vel_accum = body(3).err_vel_accum + body(3).err_vel*h;
	body(3).Tc_vel = body(3).Kp_vel*body(3).err_vel + body(3).Kd_vel*(body(3).err_vel - body(3).err_vel_prev)/h + body(3).Ki_vel*body(3).err_vel_accum;
	body(3).err_vel_prev = body(3).err_vel;

	body(3).T_control = body(3).Tc_pos + body(3).Tc_vel;
    body(3).T_control = 0;

	body(3).Ki = body(3).Mih;
	body(2).Ki = body(2).Mih + body(3).Ki;
	body(1).Ki = body(1).Mih + body(2).Ki;
	
	body(1).Di_temp = body(1).Di;
	body(2).Di_temp = body(1).Di + body(2).Di;
	body(3).Di_temp = body(1).Di + body(2).Di + body(3).Di;

 	body(1).Tg = body(1).Bi'*(body(1).Fg - body(1).Ki*body(1).Di_temp);
 	body(2).Tg = body(2).Bi'*(body(2).Fg - body(2).Ki*body(2).Di_temp);
 	body(3).Tg = body(3).Bi'*(body(3).Fg - body(3).Ki*body(3).Di_temp);
    
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
        
        body(i).yp = body(i).Ta - body(i).r_hat;
    end

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