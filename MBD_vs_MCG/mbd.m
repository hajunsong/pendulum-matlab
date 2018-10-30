clc; clear all; close all;
global  AT AT1
format long g

t_end = 5;
t = 0;
h = 0.001;

g = -9.80665;

%% read_basebody
A0 = [0,0,1;1,0,0;0,1,0];
C01 = eye(3);
s01p = [0;0;0];

%% read_body
body1.qi = 0;
body1.dqi = 0;

body1.ri = [0;0;0];

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

body2.qi = 0;
body2.dqi = 0;

body2.ri = [0;0;0];

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

body = [body1, body2];

num_body = 2;

%% define Y vector
if num_body == 1
    Y(1,1) = body.qi;
    Y(2,1) = body.dqi;
else
    for i = 1 : num_body
        Y(i,1) = body(i).qi;
    end
    for i = 1 : num_body
        Y(i+num_body,1) = body(i).dqi;
    end
end
        
Yp = zeros(2*num_body,1);
Yp_old = Yp;
index = 1;
intcount = 1;

%% mbd
while(t <= t_end)
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
        
        body(i).Qih_g = [body(i).Fic; body(i).rict*body(i).Fic];
        body(i).Qih_c = [body(i).mi*body(i).drict*body(i).wi; body(i).mi*body(i).rict*body(i).drict*body(i).wi - body(i).wit*body(i).Jic*body(i).wi];
        disp([i;1;body(i).mi*body(i).drict*body(i).wi]');
        disp([i;2;body(i).mi*body(i).rict*body(i).drict*body(i).wi]');
        disp([i;3;body(i).wit*body(i).Jic*body(i).wi]');
        body(i).Qih_a = [zeros(3,1);body(i).Tic];

        %% Velocity Coupling Body
        body(i).drit = tilde(body(i).dri);
        if i == 1
            body(i).dHi = zeros(3,1);
        else
            body(i).dHi = body(i-1).wit*body(i).Hi;
        end
        body(i).Di = [body(i).drit*body(i).Hi + body(i).rit*body(i).dHi; body(i).dHi]*body(i).dqi;
    end
    
    %% system EQM
    for i = num_body : -1 : 1
        body(i).Ki = body(i).Mih;
        body(i).Li = body(i).Qih;
        body(i).Li_g = body(i).Qih_g;
        body(i).Li_c = body(i).Qih_c;
        body(i).Li_a = body(i).Qih_a;
        if i ~= num_body
            body(i).Ki = body(i).Ki + body(i+1).Ki;
            body(i).Li = body(i).Li + body(i+1).Li - body(i+1).Ki*body(i+1).Di;
            body(i).Li_g = body(i).Li_g + body(i+1).Li_g - body(i+1).Ki*body(i+1).Di;
            body(i).Li_c = body(i).Li_c + body(i+1).Li_c - body(i+1).Ki*body(i+1).Di;
            body(i).Li_a = body(i).Li_a + body(i+1).Li_a - body(i+1).Ki*body(i+1).Di;
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
    end
    
    for i = 1 : num_body
        D_temp = zeros(6,1);
        for j = 1 : i
            D_temp = D_temp + body(j).Di;
        end

        body(i).Tg = body(i).Bi'*(body(i).Li_g - body(i).Ki*D_temp);
        body(i).Tc = body(i).Bi'*(body(i).Li_c - body(i).Ki*D_temp);
        body(i).Ta = body(i).Bi'*(body(i).Li_a - body(i).Ki*D_temp);
  
        Q(i,1) = body(i).Tg + body(i).Tc + body(i).Ta;
%         disp([body(i).Tg, body(i).Tc, body(i).Ta])
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
    
    %% integration
%     Y = Y + Yp_old*h + 0.5*h^2*(Yp - Yp_old);
    [Y, t_next, intcount] = absh3(t, Y, Yp, h, intcount);
        
    data(index,1:3) = [t, body(1).qi, body(2).qi];
%     disp(data(index,:));
    
    Yp_old = Yp;
%     t = t + h;
    t = t_next;
    index = index + 1;
end

data2 = load(sprintf('body%d.txt',num_body));
figure
set(gcf,'Color',[1,1,1])
subplot(2,1,1)
plot(data(:,1), data(:,2),'LineWidth',2)
hold on
plot(data2(:,1), data2(:,2),'--','LineWidth',2)
grid on
xlabel('Time [s]')
ylabel('Position [rad]')
legend('MATLAB','ADAMS')

subplot(2,1,2)
plot(data(:,1), data(:,3),'LineWidth',2)
hold on
plot(data2(:,1), data2(:,5),'--','LineWidth',2)
grid on
xlabel('Time [s]')
ylabel('Position [rad]')
legend('MATLAB','ADAMS')