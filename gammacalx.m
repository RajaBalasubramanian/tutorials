        
               ... Sine wave generation using Extended Taylor Series algorithm
                        ...Code for calculation of gamma values 

clc;
clear all;
close all;

% Digital equivalent of 0 to 45 is 0 to 1608
% Number of values for 0-360 is 45*8 
% i.e., number of values for 0 to (2*pi) is 1608*8= 12864
% Resolution is ((2*pi)/12864)=4.8843e-004

%  t=0.000001:0.0175:0.7854;
t=0.000001:4.8843e-004:0.7854; % pi/4= 0.7854 ; 4.8843e-004 is resolution

                             ... Gamma calculation for sine

% Proposed sine computation polynomial is: sin(t)=(t-((u.^2/2).*t))
% u is a multiplicative constant in proposed sine computation polynomial 
% u(t)=sqrt(2*(t-sin(t))./t)

ui=sqrt(2*(t-sin(t))./t); 
% ui refers to ideal gamma values for sine determined using u(t)=sqrt(2*(t-sin(t))./t)

%u1 refers to fixed increment gamma values for sine
% u1 corresponding to t=0.7854 is chosen ; u(0.7854)=0.4465
% u1 is incremented uniformly over the range 0-1608 with step size of 0.4465/1608
u1=0:0.4465/1608:0.4465; 

figure,plot(t,ui,'r',t,u1,'b'); 
% Plot for deviation in ideal and fixed increment gamma values for sine 
legend('ideal','practical');
ylabel('Gamma for sine');
xlabel('Phase value (radians)');
title('Constant for sine');

% errors refers to deviation in gamma values for sine 
errors=ui-u1;  % Difference between ideal and fixed increment gamma values for sine 
figure,plot(t,errors);
% Plot for difference between ideal and fixed increment gamma values for cosine
title('Error for sine constant'); 

                            ... Gamma calculation for cosine
                                
% Proposed cosine computation polynomial is: cos(t)=(1-(t.*r)+((r.*r)/2))
% r is a multiplicative constant in proposed cosine computation polynomial
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

figure,plot(t,ri,'r',t,r1,'b');
% Plot for deviation in ideal and fixed increment gamma values for cosine 
legend('ideal','practical');
ylabel('Gamma for cosine');
xlabel('Phase value (radians)');
title('Constant for cosine');

% errorc refers to deviation in gamma values for cosine 
errorc=ri-r1; % Difference between ideal and fixed increment gamma values for cosine
figure,plot(t,errorc); 
% Plot for difference between ideal and fixed increment gamma values for cosine
title('Error for cosine constant'); 

