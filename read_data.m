function read_data

    global start_time end_time h g
    global base body1 body2 body3 body4 body5 body6
    global r_hat r_hat_dot K A0 C01 s01p
    global des_pos err_pos err_pos_accum err_pos_prev Kp_pos Ki_pos Kd_pos
    global des_vel err_vel err_vel_accum err_vel_prev Kp_vel Ki_vel Kd_vel
    
    %% read_system
    start_time = 0;
    end_time = 5;
    h = 0.001;
    
    g = -9.80665;
    
    %% read_basebody
    
    A0 = [0,0,1;1,0,0;0,1,0];
    C01 = eye(3);
    s01p = [0;0;0];
    
    %% read_body1
    body1.qi = 0;
    body1.dqi = 0;
    
    body1.ri = [0;0;0];
    body1.pi = [0;0;0;0];
    
    body1.dri = [0;0;0];
    body1.wi = [0;0;0];
    
    body1.mi = 2;
    
    body1.Jip = [1.5 0 0;
        0 1.5 0;
        0 0 1.5];
    
    body1.rhoip = [0.15; 0; 0];
    
    body1.Cii = [0 0 -1;
        0 1 0;
        1 0 0];
    
    body1.Cij = [0 -1 0;
        1 0 0;
        0 0 1];
    
    body1.sijp = [0.3;0;0];
    
    body1.r_hat = 0;
    body1.r_hat_dot = 0;
    body1.K = 500;
    
    %% read_body2
    body2.qi = 0;
    body2.dqi = 0;
    
    body2.ri = [0;0;0];
    body2.pi = [0;0;0;0];
    
    body2.dri = [0;0;0];
    body2.wi = [0;0;0];
    
    body2.mi = 5;
    
    body2.Jip = [3 0 0;
        0 3 0;
        0 0 3];
    
    body2.rhoip = [0; -0.2; 0];
    
    body2.Cii = [0 1 0;
        0 0 1;
        1 0 0];
    
    body2.Cij = [0 -1 0;
        1 0 0;
        0 0 1];
    
    body2.sijp = [0;-0.4;0];
    
    body2.r_hat = 0;
    body2.r_hat_dot = 0;
    body2.K = 500;
    
    %% read_body3
    body3.qi = 0;
    body3.dqi = 0;
    
    body3.ri = [0;0;0];
    body3.pi = [0;0;0;0];
    
    body3.dri = [0;0;0];
    body3.wi = [0;0;0];
    
    body3.mi = 10;
    
    body3.Jip = [5 0 0;
        0 5 0;
        0 0 5];
    
    body3.rhoip = [-0.25; 0; 0];
    
    body3.Cii = [0 0 1;
        0 -1 0;
        0 0 1];
    
    body3.Cij = [0 -1 0;
        1 0 0;
        0 0 1];
    
    body3.sijp = [-0.5;0;0];
    
    body3.r_hat = 0;
    body3.r_hat_dot = 0;
    body3.K = 500;
    
    %% read_body4
    body4.qi = 0;
    body4.dqi = 0;
    
    body4.ri = [0;0;0];
    body4.pi = [0;0;0;0];
    
    body4.dri = [0;0;0];
    body4.wi = [0;0;0];
    
    body4.mi = 15;
    
    body4.Jip = [7 0 0;
        0 7 0;
        0 0 7];
    
    body4.rhoip = [0;0.3;0];
    
    body4.Cii = [0 -1 0;
        0 0 -1;
        1 0 0];
    
    body4.Cij = [0 -1 0;
        1 0 0;
        0 0 1];
    
    body4.sijp = [0;0.6;0];
    
    body4.r_hat = 0;
    body4.r_hat_dot = 0;
    body4.K = 500;
    
    %% read_body5
    body5.qi = 0;
    body5.dqi = 0;
    
    body5.ri = [0;0;0];
    body5.pi = [0;0;0;0];
    
    body5.dri = [0;0;0];
    body5.wi = [0;0;0];
    
    body5.mi = 20;
    
    body5.Jip = [10 0 0;
        0 10 0;
        0 0 10];
    
    body5.rhoip = [0.35;0;0];
    
    body5.Cii = [0 0 -1;
        0 1 0;
        1 0 0];
    
    body5.Cij = [0 -1 0;
        1 0 0;
        0 0 1];
    
    body5.sijp = [0.7;0;0];
    
    body5.r_hat = 0;
    body5.r_hat_dot = 0;
    body5.K = 500;
    
    %% read_body6
    body6.qi = 0;
    body6.dqi = 0;
    
    body6.ri = [0;0;0];
    body6.pi = [0;0;0;0];
    
    body6.dri = [0;0;0];
    body6.wi = [0;0;0];
    
    body6.mi = 25;
    
    body6.Jip = [12 0 0;
        0 12 0;
        0 0 12];
    
    body6.rhoip = [0;-0.4;0];
    
    body6.Cii = [0 1 0;
        0 0 1;
        1 0 0];
    
    body6.Cij = [0 -1 0;
        1 0 0;
        0 0 1];
    
    body6.sijp = [0;-0.8;0];
    
    body6.r_hat = 0;
    body6.r_hat_dot = 0;
    body6.K = 500;
    
    %% read_control
    des_pos = pi/180;
    err_pos = 0;
    err_pos_accum = 0;
    err_pos_prev = 0;
    
    Kp_pos = 30000;
    Ki_pos = 0;
    Kd_pos = 2;
    
    des_vel = 0;
    err_vel = 0;
    err_vel_accum = 0;
    err_vel_prev = 0;
    
    Kp_vel = 100;
    Ki_vel = 0;
    Kd_vel = 0;
    
    body1.des_pos = pi;
    body1.err_pos = 0;
    body1.err_pos_accum = 0;
    body1.err_pos_prev = 0;
    
    body1.Kp_pos = 30000;
    body1.Ki_pos = 0;
    body1.Kd_pos = 2;
    
    body1.des_vel = 0;
    body1.err_vel = 0;
    body1.err_vel_accum = 0;
    body1.err_vel_prev = 0;
    
    body1.Kp_vel = 100;
    body1.Ki_vel = 0;
    body1.Kd_vel = 0;
    
    body2.des_pos = 0*pi/180;
    body2.err_pos = 0;
    body2.err_pos_accum = 0;
    body2.err_pos_prev = 0;
    
    body2.Kp_pos = 30000;
    body2.Ki_pos = 0;
    body2.Kd_pos = 2;
    
    body2.des_vel = 0;
    body2.err_vel = 0;
    body2.err_vel_accum = 0;
    body2.err_vel_prev = 0;
    
    body2.Kp_vel = 100;
    body2.Ki_vel = 0;
    body2.Kd_vel = 0;
        
    body3.des_pos = 0*pi/180;
    body3.err_pos = 0;
    body3.err_pos_accum = 0;
    body3.err_pos_prev = 0;
    
    body3.Kp_pos = 30000;
    body3.Ki_pos = 0;
    body3.Kd_pos = 2;
    
    body3.des_vel = 0;
    body3.err_vel = 0;
    body3.err_vel_accum = 0;
    body3.err_vel_prev = 0;
    
    body3.Kp_vel = 100;
    body3.Ki_vel = 0;
    body3.Kd_vel = 0;
    
    body4.des_pos = 0*pi/180;
    body4.err_pos = 0;
    body4.err_pos_accum = 0;
    body4.err_pos_prev = 0;
    
    body4.Kp_pos = 30000;
    body4.Ki_pos = 0;
    body4.Kd_pos = 2;
    
    body4.des_vel = 0;
    body4.err_vel = 0;
    body4.err_vel_accum = 0;
    body4.err_vel_prev = 0;
    
    body4.Kp_vel = 100;
    body4.Ki_vel = 0;
    body4.Kd_vel = 0;
    
    body5.des_pos = 0*pi/180;
    body5.err_pos = 0;
    body5.err_pos_accum = 0;
    body5.err_pos_prev = 0;
    
    body5.Kp_pos = 30000;
    body5.Ki_pos = 0;
    body5.Kd_pos = 2;
    
    body5.des_vel = 0;
    body5.err_vel = 0;
    body5.err_vel_accum = 0;
    body5.err_vel_prev = 0;
    
    body5.Kp_vel = 100;
    body5.Ki_vel = 0;
    body5.Kd_vel = 0;
    
    body6.des_pos = 0*pi/180;
    body6.err_pos = 0;
    body6.err_pos_accum = 0;
    body6.err_pos_prev = 0;
    
    body6.Kp_pos = 30000;
    body6.Ki_pos = 0;
    body6.Kd_pos = 2;
    
    body6.des_vel = 0;
    body6.err_vel = 0;
    body6.err_vel_accum = 0;
    body6.err_vel_prev = 0;
    
    body6.Kp_vel = 100;
    body6.Ki_vel = 0;
    body6.Kd_vel = 0;
    
end