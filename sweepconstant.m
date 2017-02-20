                 
                   ... Sine wave generation for FCW=1 
                     
...Sweep constant to minimise hardware overhead (ROM size)
... Constant values vary uniformly for each increment of phase accumulator output 
... Angle values are in terms of radians

clc;
clear all;
close all;

% Digital equivalent of 0 to 45 is 0 to 1608
% Number of values for 0-360 is 45*8 
% i.e., number of values for 0 to (2*pi) is 1608*8= 12864
% Resolution is ((2*pi)/12864)=4.8843e-004

                   ...Plot for proposed sin(0:pi/4)

t=0.000001:4.8843e-004:0.7854;   % pi/4= 0.7854 ; 4.8843e-004 is resolution
% Proposed sine computation polynomial is: sin(t)=(t-((u.^2/2).*t))
% u is a multiplicative constant in proposed sine computation polynomial 
% u(t)=sqrt(2*(t-sin(t))./t);
% u corresponding to t=0.7854 is chosen ; u(0.7854)=0.4465
% u is incremented uniformly over the range 0-1608 with step size of 0.4465/1608
u=0:0.4465/1608:0.4465;
sina=(t-((u.^2/2).*t)); %sine values are calculated 
% comparative plot for ideal and proposed sine(0 to pi/4)
figure,plot(t,sina,'b',t,sin(t),'black'); 
legend('proposed','ideal');
xlabel('Angle (radians)');
ylabel(' Amplitude ');
title('Comparison Plot for sine ( 0 : pi/4 )');

                    ...Plot for proposed cos(0:pi/4)

t=0.000001:4.8843e-004:0.7854;   % pi/4= 0.7854 ; 4.8843e-004 is resolution
% Proposed cosine computation polynomial is: cos(t)=(1-(t.*r)+((r.*r)/2))
% r is a multiplicative constant in proposed cosine computation polynomial
% r corresponding to t=0.7854 is chosen ; r(0.7854)=0.6091
% r is incremented uniformly over the range 0-1608 with step size of 0.6091/1608 
r=0:0.6091/1608:0.6091; 
cosa=(1-(t.*r)+((r.*r)/2)); % cosine values are determined
% comparative plot for ideal and proposed cosine(0 to pi/4)
figure,plot(t,cosa,'r',t,cos(t),'black');
legend('proposed','ideal');
xlabel('Angle (radians)');
ylabel(' Amplitude ');
title('Comparison Plot for cos ( 0 : pi/4 )');

                  ...Plot for proposed sin(0:(2*pi))
                     ... Sine-Cosine symmetry is used
                         ... Frequency Control word=1
out=[sina fliplr(cosa(1:1608)) cosa(2:1609) fliplr(sina(1:1608)) -sina(2:1609) -fliplr(cosa(1:1608)) -cosa(2:1609) -fliplr(sina(1:1608)) ];
l=0:(pi*2/12864):pi*2;
% comparative plot for ideal and proposed sin(0 to(2*pi))
figure,plot(l,out,'r',l,sin(l),'b');
legend('proposed','ideal');
xlabel('Angle (radians)');
ylabel(' Amplitude ');
title('Plot for sine ( 0 : (2*pi) )');

fin = 1;
 fs = 12864;
  N = 12864;
%  fs=12865;
%  N=12865;
adc_resolution = 12; % Or Quantizer resolution


Nw = floor((fin*N)/fs);
amp = 0.449;
%  signal = amp*sin(2*pi*Nw/N*[0:N-1]);
 
signal=[sina fliplr(cosa(1:1608)) cosa(2:1609) fliplr(sina(1:1608)) -sina(2:1609) -fliplr(cosa(1:1608)) -cosa(2:1609) -fliplr(sina(1:1608)) ];

figure(1);
subplot (2,1,1);
plot(signal);
title('Original Signal');
ylabel('Signal Amplitude');
xlabel('Time');

vref=2;
res=vref/(power(2,adc_resolution)-1);
quantize_signal=round(signal/res);
% Scalar quantization is implemented over here
% no_quantiz_levels = power(2, adc_resolution);
% quantize_signal=floor((no_quantiz_levels-1)*signal)/(no_quantiz_levels/2);

subplot (2,1,2);
plot(signal-quantize_signal/2, 'g*');
title('Quantization error');

% FFT and SNR measurement
% Remember that hann window reduces the amplitude of the signal to half in
% frequency domain and that is why scaling is done by N/4 instead of N/2
window_type=ones(1,N);%hann(N)'; %Coherent sampling, no need for windowing
window_scaling_factor=sum(window_type);
% fft_signal = fft(quantize_signal.*window_type,N)*(2/window_scaling_factor);

fft_signal=fft(quantize_signal,N);

freq_scal = linspace (0, 0.5, N/2); % for normalized frequency
X=fft_signal.*conj(fft_signal);

%  signal_value=(sum(fft_signal(1:N/2).*conj(fft_signal(1:N/2))))
% signal_values = 2*(signal_value)
figure(2);
plot(freq_scal, 10*log10(X(1:N/2)), '-x');
hold on
% plot(freq_scal(Nw:Nw+2), 10*log10(abs(fft_signal(Nw:Nw+2))), 'r*');
plot( 10*log10(abs(fft_signal(Nw:Nw+2))), 'r*');
axis([0 0.01 -120 10]);
title('FFT Plot (N = 12864, Hann window is used)');
xlabel('Frequency (Normalized)');

signal_value = (sum(fft_signal(Nw:Nw+2).*conj(fft_signal(Nw:Nw+2))))

noise_bins = [fft_signal(2:fin-1) fft_signal(fin+3:N/2)];
 


noise = (sum((abs(noise_bins)).^2)) /((N/2)-3)
noise = 2*(sum(noise_bins.*conj(noise_bins)))

noise_formula = ((1/power(2,adc_resolution))^2)/12

% SNR = 10*log10(signal_value/noise)
% SNR_byQuantizationFormula = 10*log10(signal_value/noise_formula)

                  ... Spectrum plot

N=12864;  % number of values is 12864
fs=12864; 
f=fs*(0:(N-1)/2)/N; % freuency bin f=fs/N 
sp=abs(fft(out,N)); % Taking fft and considering absolute values
% signal_value=sum(abs(fft(out,N)))
sp=abs(fft(quantize_signal,N));
yfft=(sp)/max(sp);  % Normalise fft vaues with amplitude of fundamental bin
figure,plot(f,20*log(yfft(1:(N/2)))); % Plot of power spectrum in dB
axis([0,100,-500,0]); % Observe a specified range of spectrum
xlabel('Frequency bins (Hz)');
ylabel('Power (dB)');
title(' Spectrum for proposed sine computation');

                  ...To display sfdr value
%sfdr is range between fundamental and highest spur
sfdrv=sort(abs(20*log(yfft(1:N/2)))); % sort sfdr values in ascending order
sfdrv=sfdrv(2:end); % Discard fundamental bin's power value i.e. 0
sfdr=min(sfdrv); 
disp('sfdr (in dBc) is:');
disp(sfdr); % Display sfdr value
