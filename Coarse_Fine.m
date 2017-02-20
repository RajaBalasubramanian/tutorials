clc;
close all;
clear all;

coarse=12;
fine=3;
N=coarse+fine;
dac=N+1;
vref = 2;

dac_res=vref/(2^dac-1);

a1=floor((2*pi)*2^coarse);
b1=floor(a1/8);
a1=b1*8;
res1 = (pi/4)/b1;

a2=floor((2*pi)*2^fine);
b2=floor(a2/8);
a2=b2*8;
res2 = res1/b2;

l1=0:res1:pi/4;
sina1=sin(l1);
cosa1=cos(l1);
 
l2=0:res2:res1;
sina2=sin(l2);
cosa2=cos(l2);

b=((b1+1)*(b2+1))-1;
a=b*8;
res=(pi/4)/b;

l=0:res:(pi/4);

k=1;
for i=1:1:b1+1
    for j=1:1:b2+1
        sina(k)=(sina1(i)*cosa2(j))+(cosa1(i)*sina2(j));
        cosa(k)=(cosa1(i)*cosa2(j))-(sina1(i)*sina2(j));
        k=k+1;
    end
end

t=0:res:(2*pi);
signal=[sina fliplr(cosa(1:b)) cosa(2:b+1) fliplr(sina(1:b)) -sina(2:b+1) -fliplr(cosa(1:b)) -cosa(2:b+1) -fliplr(sina(1:b)) ];

figure,plot(t,signal); 
xlabel('Angle (radians)');
ylabel(' Amplitude ');
title('Plot for sine ( 0 : 2*pi )');

fs = a;
samples = a;

fin=1;
Nw=floor(fin*samples)/fs;
 
quantize_signal=dac_res*(round(signal/dac_res));

figure,plot(t,quantize_signal); 
xlabel('Angle (radians)');
ylabel(' Amplitude ');
title('Plot for Quantized sine ( 0 : 2*pi )');

error=signal-quantize_signal;
figure, plot(t,error);
title('Quantization error');

original = sin(t);
figure; plot(t,original,'r',t,quantize_signal,'black');
legend('ideal','proposed');
xlabel('Angle (radians)');
ylabel(' Amplitude ');
title('Comparison Plot for sine ( 0 : 2*pi )');

error = original-quantize_signal;
figure; plot(t,error);
xlabel('Angle (radians)');
ylabel(' Amplitude ');
title('Error');

window_type=ones(1,samples);
window_scaling_factor=sum(window_type);

fft_signal=fft(quantize_signal,samples);

freq_scal = linspace (0, 0.5, samples/2);
X=fft_signal.*conj(fft_signal);

figure;
plot(freq_scal, 10*log10(X(1:samples/2)), '-x');
hold on
plot( 10*log10(abs(fft_signal(Nw:Nw+2))), 'r*');
axis([0 0.01 -120 10]);
title('FFT Plot (N = 12864, Hann window is used)');
xlabel('Frequency (Normalized)');

signal_value = (sum(fft_signal(Nw:Nw+2).*conj(fft_signal(Nw:Nw+2))));
noise_bins = [fft_signal(2:fin-1) fft_signal(fin+3:samples/2)]; 
noise = (sum((abs(noise_bins)).^2)) /((samples/2)-3);
noise_value = 2*(sum(noise_bins.*conj(noise_bins)));

noise_formula = ((1/power(2,dac))^2)/12;

f=fs*(0:(samples-1)/2)/samples; 
sp=abs(fft(quantize_signal,samples));
yfft=(sp)/max(sp); 
figure,plot(f,20*log(yfft(1:(samples/2)))); 
axis([0,100,-500,0]);
xlabel('Frequency bins (Hz)');
ylabel('Power (dB)');
title(' Spectrum for proposed sine computation');

sfdrv=sort(abs(20*log(yfft(1:samples/2)))); 
sfdrv=sfdrv(2:end);
sfdr=min(sfdrv); 
disp('sfdr (in dBc) is:');
disp(sfdr); 

snr_ideal = (6.02 * dac) + 1.76;
disp('Maximum Achievable SNR (in dB) is :');
disp(snr_ideal);

snr=10*log10(signal_value/noise_value);
disp('snr (in dB) is:');
disp(snr);