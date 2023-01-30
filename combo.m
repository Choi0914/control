clc, clear
close all
%% Set Parameters
% Set Simulation Time
    delta_t = 0.001;
    end_time = 5;

    sim_time =(0:0.001:5);
% Set Sine Wave
    sine_mag1 = 2.0; sine_freq1 = 1.0; % Main Signal
    sine_mag2 = 0.5; sine_freq2 = 10.0; % Noise
%% Set FFT
    Fs = 1/delta_t; % Sampling Frequncy
    T = delta_t; % Sampling Period
    L = length(sim_time); % Length of Signal
    T_vector = (0:L-1)*T; % Time Vector

    fft_f = Fs*(0:(L/2))/L; % Frequency Range
%% define sin wave
sim_y = sine_mag1*sin(sine_freq1*(2*pi*T_vector)) + sine_mag2 * sin(sine_freq2 * (2*pi*T_vector)) + 0.8 * randn(size(T_vector));
%% MAF filter
   m1 = 1;
   n = 100;
   sum = 0;
   for t=0:delta_t:end_time
        if m1 > n
            MAF_R(m1) = MAF_R(m1-1) + (sim_y(m1)-sim_y(m1-n))/n; % 이동평균 재귀식
        else 
            sum = sum + sim_y(m1); %100개씩 묶기로 해서 MAF_R(100)까지는 재귀식을 쓸수없어 점점 더해지는 값들을 n으로 나누는 식으로 했다.
            MAF_R(m1) = sum/n;
        end
        m1 = m1+1;
   end
   %% LPF filter
   m2=2;
   tau = 0.0159;
   LPF_R(1) = sim_y(1);
   for t=delta_t:delta_t:end_time
        alpha = tau/(delta_t + tau);
        LPF_R(m2) = (1-alpha)*sim_y(m2) + alpha*LPF_R(m2-1);
        m2 = m2 + 1;
   end
%% Calc NONE FFT    
    fft_y_temp1 = abs(fft(sim_y)/L);
    fft_y1 = fft_y_temp1(1:L/2+1);
    fft_y1(2:end-1) = 2*fft_y1(2:end-1);
%% Calc MAF FFT
    fft_y_temp2 = abs(fft(MAF_R)/L);
    fft_y2 = fft_y_temp2(1:L/2+1);
    fft_y2(2:end-1) = 2*fft_y2(2:end-1);
%% Calc LPF FFT
    fft_y_temp3 = abs(fft(LPF_R)/L);
    fft_y3 = fft_y_temp3(1:L/2+1);
    fft_y3(2:end-1) = 2*fft_y3(2:end-1);
%% Draw Graph
figure('units', 'pixels', 'pos',[500 500 500 1000],'Color',[1,1,1]);
    subplot(4,1,1) % Time-Domain
        Xmin = 0.0; XTick = 1.0; Xmax = end_time;
        Ymin = -5.0; YTick = 1.0; Ymax = 5.0;

            plot(T_vector, sim_y,'-k','LineWidth',1) %None
            hold on;
            plot(T_vector, MAF_R,'-y','LineWidth',1) %MAF
            hold on;
            plot(T_vector, LPF_R,'-r','LineWidth',1) %LPF

            grid on; % Grid On
            axis([Xmin Xmax Ymin Ymax]) % Graph ## ## ##
            set(gca,'XTick',(Xmin:XTick:Xmax)); % X# Grid ##
            set(gca,'YTick',(Ymin:YTick:Ymax)); % Y# Grid ##
        xlabel('Time (s)', 'fontsize',20);
        ylabel('Magnitude', 'fontsize',20);
        title ('Time Domain', 'fontsize',25);
    
    subplot(4,1,2) % Frequncy-Domain
        Xmin = 0.0; Xmax = 11.0;
        Ymin = 0.0; Ymax = 3.0;

            stem(fft_f, fft_y1,'-k','LineWidth',2)

            grid on; % Grid On
            axis([Xmin Xmax Ymin Ymax]) % Graph ## ## ##
            set(gca,'XTick',(0.0:1.0:11.0)); % X# Grid ##
            set(gca,'YTick',(0.0:0.5:3.0)); % Y# Grid ##
        xlabel('Frequency (Hz)', 'fontsize',20);
        ylabel('Magnitude', 'fontsize',20);
        title ('None Filter', 'fontsize',25);

    subplot(4,1,3) % Frequncy-Domain
        Xmin = 0.0; Xmax = 11.0;
        Ymin = 0.0; Ymax = 3.0;

            stem(fft_f, fft_y2,'-k','LineWidth',2)

            grid on; % Grid On
            axis([Xmin Xmax Ymin Ymax]) % Graph ## ## ##
            set(gca,'XTick',[0:1.0: 11.0]); % X# Grid ##
            set(gca,'YTick',[0.0:0.5:3.0]); % Y# Grid ##
        xlabel('Frequency (Hz)', 'fontsize',20);
        ylabel('Magnitude', 'fontsize',20);
        title ('Moving Average Filter', 'fontsize',25);

    subplot(4,1,4) % Frequncy-Domain
        Xmin = 0.0; Xmax = 11.0;
        Ymin = 0.0; Ymax = 3.0;

            stem(fft_f, fft_y3,'-k','LineWidth',2)

            grid on; % Grid On
            axis([Xmin Xmax Ymin Ymax]) % Graph ## ## ##
            set(gca,'XTick',[0:1.0: 11.0]); % X# Grid ##s
            set(gca,'YTick',[0.0:0.5:3.0]); % Y# Grid ##
        xlabel('Frequency (Hz)', 'fontsize',20);
        ylabel('Magnitude', 'fontsize',20);
        title ('Low Pass Filter', 'fontsize',25);