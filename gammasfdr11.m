clc;
clear all;
close all;

%  x=diff(u1,t);

% ui=sqrt(2*(t-sin(t))./t); 
% x=diff(ui,t);

% % for i=0:1608;
% %     for j=0:1608;
% %      u(i)=0:0.4465/1608:0.4465;
% %      t(j)=0.000001:4.8843e-004:0.7854;
% %      y=(u(i)-u(i-1))/(t(j)-t(j-1));
% %     end
% % end



% t=0.000001:4.8843e-004:0.7854;
%  x=0:0.4465/1608:0.4465;
% dx = x(2:end) - x(1:end-1);
% dt = t(2:end) - t(1:end-1);
% v = dx./dt
% plot(v,t);

%  v=diff(sqrt(2*(1-(sin(t)/t))));

% % for i=1:1:1608;
% %     dx(i)=x(i+1)-x(i);
% %     dt(i)=t(i+1)-t(i);
% % end
% % v=dx./dt;
% % v1=[v,0];
% % plot(t,v1);
% % ylabel('derivative');
% % xlabel('Phase value (radians)');



a=0:0.2416/535:0.2416;
% disp(r1);

m1=0:0.2613/535:0.2613;

for i=1:1:535;
    da(i)=a(i+1)-a(i);
    dm1(i)=m1(i+1)-m1(i);
end
v1=da./dm1;
%  v11=[0,v1];
% figure,plot(a,v11);
% figure,plot(m1,a);

b=0.2420:0.2025/535:0.4445;

m2=0.2618:0.2613/535:0.5231;

for i=1:1:535;
    db(i)=b(i+1)-b(i);
    dm2(i)=m2(i+1)-m2(i);
end
v2=db./dm2;
% v22=[0,v2];
% figure,plot(b,v22);

% figure,plot(m2,b);
c=0.4448:0.1643/536:0.6091;

m3=0.5236:0.2618/536:0.7854;
r=[a b c];

for i=1:1:535;
    dc(i)=c(i+1)-c(i);
    dm3(i)=m3(i+1)-m3(i);
end
v3=dc./dm3;
% v33=[0,v3];
% figure,plot(c,v33);

% figure,plot(m3,c);

 figure,plot(m1,a,'b',m2,b,'g',m3,c,'b');
legend('linearization');
xlabel('phasevalue(radians)');
ylabel('gamma for cosine');
title('constant for cosine');

t=0.000001:4.8843e-004:0.7854;
i=1;
for t1=0.000001:4.8843e-004:0.7854;
    p=[0.5 (-t1) (1-cos(t1))]; % coefficients of the quadratic equation
ri(i)=min(roots(p)); % choose minimum of two roots for each phase value
% ri refers to ideal gamma values for cosine 
i=i+1;
end

% r1 refers to fixed increment gamma values for cosine
% r1 corresponding to t=0.7854 is chosen ; r(0.7854)=0.6091
% r1 is incremented uniformly over the range 0-1608 with step size of 0.6091/1608
r1=0:0.6091/1608:0.6091; 

figure,plot(t,ri,'r');
% Plot for deviation in ideal and fixed increment gamma values for cosine 
legend('ideal');
ylabel('Gamma for cosine');
xlabel('Phase value (radians)');
title('Constant for cosine');



figure,plot(t,ri,'r',m1,a,'b',m2,b,'g',m3,c,'b');
legend('ideal','linearization');
xlabel('phasevalue(radians)');
ylabel('gamma for cosine');
title('constant for cosine');

figure,plot(t,r1,'g');
legend('practical');
xlabel('phasevalue(radians)');
ylabel('gamma for cosine');
title('constant for cosine');



u=0:0.4465/1608:0.4465;
sina=(t-((u.^2/2).*t)); %sine values are calculated 
% comparative plot for ideal and proposed sine(0 to pi/4)
figure,plot(t,sina,'b',t,sin(t),'black'); 
legend('proposed','ideal');
xlabel('Angle (radians)');
ylabel(' Amplitude ');
title('Comparison Plot for sine ( 0 : pi/4 )');


cosa=(1-(t.*r)+((r.*r)/2)); % cosine values are determined
% % comparative plot for ideal and proposed cosine(0 to pi/4)
figure,plot(t,cosa,'r');
legend('proposed','ideal');
xlabel('Angle (radians)');
ylabel(' Amplitude ');
title('plot for cosine');



out=[sina fliplr(cosa(1:1608)) cosa(2:1609) fliplr(sina(1:1608)) -sina(2:1609) -fliplr(cosa(1:1608)) -cosa(2:1609) -fliplr(sina(1:1608)) ];
l=0:(pi*2/12864):pi*2;
% comparative plot for ideal and proposed sin(0 to(2*pi))
figure,plot(l,out,'r',l,sin(l),'b');
legend('proposed','ideal');
xlabel('Angle (radians)');
ylabel(' Amplitude ');
title('Plot for sine ( 0 : (2*pi) )');


N=12864;  % number of values is 12864
fs=12864; 
f=fs*(0:(N-1)/2)/N; % freuency bin f=fs/N 
sp=abs(fft(out,N)); % Taking fft and considering absolute values
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




% s=0:0.6091/1608:0.6091; 
% cosay=(1-(t.*s)+((s.*s)/2)); % cosine values are determined
% % comparative plot for ideal and proposed cosine(0 to pi/4)
% figure,plot(t,cosa,'r',t,cosay,'b');
% legend('proposed','ideal');
% title('comparison plot for cosine');
