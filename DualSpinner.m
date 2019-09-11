%% DUAL-SPIN STABILIZATION
% This program numerically simulates a spacecraft or a rigid body having the
% specified inertia matrix and initial conditions and hence plots the
% response of its angular rates and attitude (3-2-1 euler angles) in order
% to study its stability using a flywheel.

%% Spacecraft Specifications
I = [350,0,0;0,300,0;0,0,400];  % Inertia matrix of whole spacecraft
I_w = 10;   % Moment of inertia of the flywheel
w_e = (2*pi/60)*[60,0,0]';  % equilibrium angular rates
w0 = w_e + (2*pi/60)*[1,1,1]';  % initial angular velicity with disturbance
e0 = (pi/180)*[0,0,0]';   % initial attitude (3-2-1)
omega = (2*pi/60)*[350,0,0]'; % angular rate of flywheel

%% Numerical Simulation
dt = 0.001; % time-step
t_total = 300;  % total time

t = 0:dt:t_total;   % discrete time instances

w = zeros(3,length(t)); % angular rates
w(:,1) = w0;
e = zeros(3,length(t)); % euler angles
e(:,1) = e0;

for i = 1:(length(t)-1)
    % angular rates --> euler angle rates
    e_dot = [cos(e(2,i)),sin(e(1,i))*sin(e(2,i)),cos(e(1,i))*sin(e(2,i));
            0,cos(e(1,i))*cos(e(2,i)),-sin(e(1,i))*cos(e(2,i));
            0,sin(e(1,i)),cos(e(1,i))]*w(:,i)/cos(e(2,i));
    e(:,i+1) = e(:,i) + e_dot*dt;
    
    % -90 < theta < 90
    if e(2,i+1) > pi/2
        e(2,i+1) = pi - e(2,i+1);
        e(1,i+1) = pi + e(1,i+1);
        e(3,i+1) = pi + e(3,i+1);
    elseif e(2,i+1) < -pi/2
        e(2,i+1) = -pi - e(2,i+1);
        e(1,i+1) = pi + e(1,i+1);
        e(3,i+1) = pi + e(3,i+1);
    end
    
    % -180 < (phi,psi) < 180
    if e(1,i+1) > pi
        e(1,i+1) = mod(e(1,i+1),-pi);
    elseif e(1,i+1) < -pi
        e(1,i+1) = mod(e(1,i+1),-pi) + pi;
    end
    if e(3,i+1) > pi
        e(3,i+1) = mod(e(3,i+1),-pi);
    elseif e(3,i+1) < -pi
        e(3,i+1) = mod(e(3,i+1),-pi) + pi;
    end
    
    % system dynamics
    w_dot = -I\cross(w(:,i),I*w(:,i) + I_w*omega);
    w(:,i+1) = w(:,i) + w_dot*dt;
end

%% Numerical Error Estimation
%
err = 60/(2*pi)*abs(w(1,:) - w_e(1)); % w.r.t. equilibrium
ref = max(err(1:6/dt)); % first peak
err = (err - ref)/ref;  % w.r.t. first peak
err = 100*(err>0).*err; % truncated negative values
max_percent_err = [max(err(1,:))]
%
figure; % error propagated over time
plot(t,err(1,:),'LineWidth',2);
grid on;
title('Numerical Error in \omega_1');
xlabel('Time (sec)');
ylabel('% Error');
%}
%% Plotting Responses
%
figure; % angular rates response
plot(t,60/(2*pi)*w(1,:),'LineWidth',2);
hold on;
grid on;
plot(t,60/(2*pi)*w(2,:),'LineWidth',2);
plot(t,60/(2*pi)*w(3,:),'LineWidth',2);
title('Angular Rates');
xlabel('Time (sec)');
ylabel('Angular Rates, \omega (rpm)');
legend({'\omega_1','\omega_2','\omega_3'},'FontSize',12);
%
figure; % attitude response
plot(t(1:(end-1)/20),180/(pi)*e(1,1:(end-1)/20),'.','LineWidth',2);
hold on;
grid on;
plot(t(1:(end-1)/20),180/(pi)*e(2,1:(end-1)/20),'.','LineWidth',2);
plot(t(1:(end-1)/20),180/(pi)*e(3,1:(end-1)/20),'.','LineWidth',2);
title('Attitude');
xlabel('Time (sec)');
ylabel('Euler Angles (deg)');
legend({'\phi','\theta','\psi'},'FontSize',12);
%}