
                 ... Sine wave generation
                     
           ... Angle values are in terms of radians
      ... Single constant value is used in sin/cos computation expression
      
clc ;
clear all;
close all;

% Input the frequency control word
in=input('enter the cw: ');
in=in*(pi/(45*4)); % Convert FCW in degrees to radians
inx=in;

% Proposed sine computation polynomial is: sin(t)=(t-((u1.^2/2).*t))
% u1 is a multiplicative constant in proposed sine computation polynomial 
u1=0.0800; 
% Proposed cosine computation polynomial is: cos(t)=(1-(t.*u2)+((u2.*u2)/2))
% u2 is a multiplicative constant in proposed cosine computation polynomial
u2=0.3500;

              ... Plot for sin( 0:(2*pi))

% Digital equivalent of 0 to 45 is 0 to 1608
% Number of values for 0-360 is 45*8 
% i.e., number of values for 0 to (2*pi) is 1608*8
i=1;
j=1;
% Sine-Cosine symmetry from 0-45 is used
for in=0:in:((1608/(2^11))*8) % radian equivalent for in=0:in:360
    if (in>0)&&(in<=(1608/(2^11))) %(in>0)&&(in<=45)
        out(j)=(in-(u1.*in));      % sin expresion is used
        j=j+1;
    end
    
    
    if (in>(1608/(2^11)))&&(in<=((1608/(2^11))*2)) %(in>45)&&(in<=90)
        out(j)=1-(u2.*(((1608/(2^11))*2)-in))+(u2.*u2/2); % cosine expresion is used
        j=j+1;
    end
    
    if (in>(1608/(2^11))*2)&&(in<=((1608/(2^11))*3))%(in>90)&&(in<=135)
        out(j)=1-(u2.*(in-((1608/(2^11))*2)))+(u2.*u2/2);  % cosine expresion is used
        j=j+1; 
    end
     
     if  (in>((1608/(2^11))*3))&&(in<=((1608/(2^11))*4))%(in>135)&&(in<=180)
        out(j)=((((1608/(2^11))*4)-in)-(u1.*(((1608/(2^11))*4)-in))); % sin expresion is used
        j=j+1;
     end
      
     if (in>((1608/(2^11))*4))&&(in<=((1608/(2^11))*5)) %(in>180)&&(in<=225)
          out(j)=(-1)*((in-((1608/(2^11))*4))-(u1.*(in-((1608/(2^11))*4)))); % sin expresion is used
        j=j+1;
     end
      
      if (in>((1608/(2^11))*5))&&(in<=((1608/(2^11))*6)) %(in>225)&&(in<=270)
         out(j)=(-1)*(1-(u2.*(((1608/(2^11))*6)-in))+(u2.*u2/2));   % cosine expresion is used
         j=j+1;
      end
      
     if (in>((1608/(2^11))*6))&&(in<=((1608/(2^11))*7)) %(in>270)&&(in<=315)
         out(j)=((-1)*(1-(u2.*(in-((1608/(2^11))*6)))+(u2.*u2/2))); % cosine expresion is used
         j=j+1;
     end
       
     if (in>((1608/(2^11))*7))&&(in<=((1608/(2^11))*8)) %(in>315)&&(in<=360)  
         out(j)=((-1)*((((1608/(2^11))*8)-in)-(u1.*(((1608/(2^11))*8)-in)))); % sin expresion is used
        j=j+1;
     end
 
end
out1=[0 out 0]; % sin values for 0 to 360
figure,plot(out1);  %Plot for sin( 0:(2*pi))
xlabel('Angle (radians)');
ylabel(' Amplitude ');
title('Plot for sine ( 0 : (2*pi) )');

% Determination of ideal sin values
isine1=sin(0:inx:((1608/(2^11))*8)); 
isin=[isine1 0];

% comparative plot for ideal and proposed sine(0 to (2*pi))
...figure,plot(0:inx:((1608/(2^11))*8),out1,'b',0:inx:((1608/(2^11))*8),isin,'black');
figure,plot(0:length(out1)-1,out1,'r',0:length(isin)-1,isin,'black');
xlabel('Angle (radians)');
ylabel(' Amplitude ');
title('Comparative Plot for sine ( 0 : (2*pi) )');

% Determination of error (deviation from ideal values)
error=isin-out1; % ideal values-calculated values
disp('Overall max. deviation:');
disp(max(error)); % Display error
figure,plot(error);
xlabel('Angle (radians)');
ylabel(' Amplitude ');
title(' Error plot');

                   ... Spectrum plot

fs=360; % sampling frequency
N=360;   % Number of values is 360 for FCW=1  
f=fs*(0:(N-1)/2)/N; % freuency bin f=fs/N 
sp=abs(fft(out,N)); % Taking fft and considering absolute values
yfft=(sp)/max(sp);  % Normalise fft vaues with amplitude of fundamental bin
figure,plot(f,20*log(yfft(1:(N/2)))); % Plot of power spectrum in dB
xlabel('Frequency bins (Hz)');
ylabel('Power (dB)');
title(' Spectrum for proposed sine computation');

                   ...To display sfdr value
%sfdr is range between fundamental and highest spur
sfdrv=sort(abs(20*log(yfft(1:N/2))));  % sort sfdr values in ascending order
sfdrv=sfdrv(2:end);% Discard fundamental bin's power value i.e. 0
sfdr=min(sfdrv);
disp('sfdr (in dBc) is:');
disp(sfdr); % Display sfdr value

pl1=out1(1:46);
pl2=out1(47:47+46);
figure,plot(0:1:45,pl1,'r',0:1:45,sind(0:45),'b');
% figure,plot(pl2,'r',,'b');

