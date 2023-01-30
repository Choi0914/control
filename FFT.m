clc, clear
close all
%%
end_time = 5;
delta_t = 0.001;
sine_mag=2;
sine_freq=1;

sine_mag2 = 0.5;
sine_freq2 = 10;
%%
Fs= 1/delta_t;
T=delta_t;
L=10000;
T_vector=(0:L-1)*T;

fft_f = Fs*(0:(L/2))/L;
%%
sin_y =sine_mag *sin(sine_freq*(2*pi*T_vector))+sine_mag2*sin(sine_freq2*(2*pi*T_vector))+0.8*randn(size(T_vector));
%%
fft_y_temp = abs(fft(sin_y)/L);
fft_y = fft_y_temp(1:L/2+1);
fft_y(2:end-1)=2*fft_y(2:end-1);

%%

figure('units', 'pixels', 'pos', [100 100 450 400], 'Color', [1,1,1]);
    subplot(2,1,1)
        Xmin=0.0;   XTick=1.0;  Xmax=end_time;
        Ymin =-3.0; YTick=1.0;  Ymax=3.0;
        plot(T_vector,sin_y, '-k', 'LineWidth', 2)

        grid on;
        axis([Xmin Xmax Ymin Ymax])
        set(gca,'XTick', [Xmin:XTick:Xmax]);
        set(gca,'YTick', [Ymin:YTick:Ymax]);
        xlabel('Time (s)', 'fontsize',20);
        ylabel('Magnitude', 'fontsize',20);
        title('Sine Wave', 'fontsize',25);
    subplot(2,1,2)
        Xmin=0.0;   Xmax=11;
        Ymin =-0.0; Ymax=3.0;
        stem(fft_f,fft_y,'-k', 'LineWidth',2)

        grid on;
        axis([Xmin Xmax Ymin Ymax])
        set(gca,'XTick', [0 1.0 10.0]);
        set(gca,'YTick', [0 0.5 2.0]);
        xlabel('Frequency (Hz)', 'fontsize',20);
        ylabel('Magnitude', 'fontsize',20);
        title('Frequency Domain', 'fontsize',25);