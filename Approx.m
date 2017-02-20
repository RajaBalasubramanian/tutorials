clc;
close all;
clear all;

N = 11;
dac = N+1;
vref = 2;

dac_res=vref/(2^dac-1);
a=floor((2*pi)*2^N);
b=floor(a/8);
a=b*8;
res = (2*pi)/a;

t=0:res:pi/4;
signal1=sin(t);
signal2=t;
error1=signal1-signal2;
signal3=1-((t.*t)/2);
signal4=cos(t);
error2=signal4-signal3;

figure, plot(t,signal1,'r',t,signal2,'g');
legend('ideal','approx');
xlabel('Angle (radians)');
ylabel(' Amplitude ');
title('Comparison Plot for sine ( 0 : pi/4 )');

figure, plot(t,signal4,'r',t,signal3,'g');
legend('ideal','approx');
xlabel('Angle (radians)');
ylabel(' Amplitude ');
title('Comparison Plot for cosine ( 0 : pi/4 )');

figure, plot(t,error1,'r',t,error2,'g');
legend('Sine error','Cosine error');
xlabel('Angle (radians)');
ylabel(' Amplitude ');
title('Approx error for sin and cos ( 0 : pi/4 )');