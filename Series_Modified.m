clc;
close all;
clear all;

N = 7;
dac = N+1;
vref = 2;

dac_res=vref/(2^dac-1);
a=floor((2*pi)*2^N);
b=floor(a/8);
a=b*8;
res = (2*pi)/a;

l=0:res:pi/4;

ui=sqrt(2*(1-(sin(l)./l)));
umax=sqrt(2*(1-(sin(pi/4)/(pi/4))));
u=0:umax/b:umax;
sina=(1-(u.^2/2)).*l;

figure,plot(l,ui,'r',l,u,'b'); 
legend('ideal','practical');
ylabel('Gamma for sine');
xlabel('Phase value (radians)');
title('Constant for sine');

errors=ui-u;
figure,plot(l,errors);
title('Error for sine constant'); 

figure,plot(l,sina,'b',l,sin(l),'black'); 
legend('proposed','ideal');
xlabel('Angle (radians)');
ylabel(' Amplitude ');
title('Comparison Plot for sine ( 0 : pi/4 )');


p=[0.5 (-pi/4) (1-cos(pi/4))];
rmax=min(roots(p)); 

i=1;
for z=0:res:pi/4
    p=[0.5 (-z) (1-cos(z))]; 
    ri(i)=min(roots(p)); 
    i=i+1;
end

points=11;

start(1)=1;
for i=1:points-1
    step(i)=floor((b+1)/points);
    finish(i)=step(i)*i; 
    start(i+1)=finish(i)+1;
end    
step(points)=(b+1)-((floor((b+1)/points))*(points-1));
finish(points)=b+1;

temp=ri; 

w=temp(start(1)):(temp(finish(1))-temp(start(1)))/(step(1)-1):temp(finish(1));
r=w;
for i=2:points
    w=temp(start(i)):(temp(finish(i))-temp(start(i)))/(step(i)-1):temp(finish(i));
    r=[r w];
end

cosa=(1-(l.*r)+((r.*r)/2));

figure,plot(l,ri,'r',l,r,'b'); 
legend('ideal','practical');
xlabel('Phase value (radians)');
ylabel('Gamma for cosine');
title('Constant for cosine');

errorc=ri-r;
figure,plot(l,errorc);
title('Error for cosine constant'); 

figure,plot(l,cosa,'b',l,cos(l),'black'); 
legend('proposed','ideal');
xlabel('Angle (radians)');
ylabel(' Amplitude ');
title('Comparison Plot for cosine ( 0 : pi/4 )');

t=0:res:2*pi;
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

figure, plot(t,signal-quantize_signal);
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