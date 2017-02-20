clc;
clear all;
close all;


a=0:0.0382/79:0.0382;
%  disp (a);

m1=0:0.0386/79:0.0386;

b=0.0386:0.0373/79:0.0759;
% disp (b);

m2=0.0391:0.0386/79:0.0777;

c=0.0764:0.0364/79:0.1128;
% disp (c);

m3=0.0781:0.0386/79:0.1167;

d=0.1133:0.0355/79:0.1488;
% disp (d);

m4=0.1172:0.0386/79:0.1558;

e=0.1492:0.0347/79:0.1839;

m5=0.1563:0.0386/79:0.1949;

f11=0.1844:0.0338/79:0.2182;

m6=0.1954:0.0386/79:0.2340;

g=0.2186:0.0329/79:0.2515;

m7=0.2344:0.0386/79:0.2730;

h=0.2520:0.032/79:0.2840;

m8=0.2735:0.0386/79:0.3121;

i11=0.2844:0.0313/79:0.3157;

m9=0.3126:0.0386/79:0.3512;

j=0.3160:0.0304/79:0.3464;

m10=0.3517:0.0386/79:0.3903;

k=0.3468:0.0295/79:0.3763;

m11=0.3907:0.0386/79:0.4293;

l=0.3767:0.0286/79:0.4053;

m12=0.4298:0.0386/79:0.4684;

m=0.4057:0.0278/79:0.4335;

m13=0.4689:0.0386/79:0.5075;

n=0.4338:0.0269/79:0.4607;

m14=0.5080:0.0386/79:0.5466;

o=0.4611:0.0261/79:0.4872;

m15=0.5470:0.0386/79:0.5856;

p1=0.4875:0.0253/79:0.5128;

m16=0.5861:0.0386/79:0.6247;

q=0.5131:0.0244/79:0.5375;

m17=0.6252:0.0386/79:0.6638;

r11=0.5378:0.0236/79:0.5614;

m18=0.6643:0.0386/79:0.7029;

s=0.5617:0.0228/79:0.5845;

m19=0.7033:0.0386/79:0.7419;

t11=0.5848:0.0243/88:0.6091;

m20=0.7424:0.043/88:0.7854;

r=[a b c d e f11 g h i11 j k l m n o p1 q r11 s t11];

t=0.000001:4.8843e-004:0.7854;
i=1;
for t1=0.000001:4.8843e-004:0.7854;
    p=[0.5 (-t1) (1-cos(t1))]; % coefficients of the quadratic equation
ri(i)=min(roots(p)); % choose minimum of two roots for each phase value
% ri refers to ideal gamma values for cosine 
i=i+1;
end


figure(1),plot(t,ri,'r',m1,a,'b',m2,b,'g',m3,c,'y',m4,d,'b',m5,e,'g',m6,f11,'b',m7,g,'g',m8,h,'b',m9,i11,'g',m10,j,'b',m11,k,'g',m12,l,'b',m13,m,'g',m14,n,'b',m15,o,'g',m16,p1,'b',m17,q,'g',m18,r11,'b',m19,s,'g',m20,t11,'b');
legend('ideal','linearization');
xlabel('phasevalue(radians)');
ylabel('gamma for cosine');
title('constant for cosine');

% figure,plot(t,ri,'r')
% 
% figure,plot(m1,a,'b',m2,b,'g',m3,c,'y',m4,d,'b',m5,e,'g',m6,f,'b',m7,g,'g',m8,h,'b',m9,i,'g',m10,j,'b');
% 
% t=[m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 m13 m14 m15 m16 m17 m18 m19 m20];
r1=0:0.6091/1608:0.6091; 

u=0:0.4465/1608:0.4465;
sina=(t-((u.^2/2).*t));

cosa=(1-(t.*r)+((r.*r)/2)); 


out=[sina fliplr(cosa(1:1608)) cosa(2:1609) fliplr(sina(1:1608)) -sina(2:1609) -fliplr(cosa(1:1608)) -cosa(2:1609) -fliplr(sina(1:1608)) ];




% Generate a simulation signal
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

figure(2);
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

 signal_value=(sum(fft_signal(1:N/2).*conj(fft_signal(1:N/2))))
% signal_values = 2*(signal_value)
figure(3);
plot(freq_scal, 10*log10(X(1:N/2)), '-x');
hold on
plot(freq_scal(Nw:Nw+2), 10*log10(abs(fft_signal(Nw:Nw+2))), 'r*');
title('FFT Plot (N = 12864, Hann window is used)');
xlabel('Frequency (Normalized)');

% signal_value = (sum(fft_signal(Nw:Nw+2).*conj(fft_signal(Nw:Nw+2))))

noise_bins = [fft_signal(2:fin-1) fft_signal(fin+3:N/2)];
% yfft = noise_bins.*conj(noise_bins); 


noise = (sum((abs(noise_bins)).^2)) /((N/2)-3)
noise = 2*(sum(noise_bins.*conj(noise_bins)))

noise_formula = ((1/power(2,adc_resolution))^2)/12

SNR = 10*log10(signal_value/noise)
SNR_byQuantizationFormula = 10*log10(signal_value/noise_formula)


N=12864;  % number of values is 12864
fs=12864; 
f=fs*(0:(N-1)/2)/N; % freuency bin f=fs/N 
% sp=abs(fft(signal,N)); % Taking fft and considering absolute values
sp=abs(fft(quantize_signal,N));
% yfft=(X)/max(X);  % Normalise fft vaues with amplitude of fundamental bin
yfft=sp/max(sp);
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






