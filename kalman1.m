clc, clear all

%%
Pos_pred = 0;           
Vel_pred = 80;          

dt = 0.01;
sim_time    = [0:0.01:5];

v_noise = 10 * randn(size(sim_time));     

Pos_pred = Pos_pred + Vel_pred * sim_time;
Vel_pred = Vel_pred + v_noise;

vel_measure = Vel_pred;

A = [1 dt; 0 1];
H = [0 1];
Q = [1 0; 0 10];
R = 500;
x = [0;80];
P = 3 * eye(2);

for(i=1:length(sim_time))
    x_hat = A*x;
    P_hat = A*P*A' +Q;

    K = P_hat*H'*inv(H*P_hat*H'+R);
    x = x_hat+K*(Vel_pred(i) - H*x_hat);
    P = P_hat - K*H*P_hat;

    Pos(i) = x(1);
    Vel(i) = x(2);
    P_gain(i) = P(1);
    K_gain(i) = K(1);
end

figure
    plot(sim_time, Pos_pred ,'r', sim_time,Pos, 'g');
    legend('pre Position', 'kalman position')
    xlabel('time(s)')
    ylabel('position')

figure
    plot(sim_time, vel_measure,'g',sim_time, Vel,'r')
    legend('Estimated','Measured')
    xlabel('time(s)')
    ylabel('velocity')