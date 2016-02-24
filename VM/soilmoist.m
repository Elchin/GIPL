function soilmoist

%hold on
m=load('19cm.txt');
nn=length(m);
%t=m(1:300,3);v=m(1:300,4);
%v(v>0.38)=NaN
figure(1)
plot(1:nn,m(:,2),'r-.')
xlabel 'time (days)'
ylabel 'volumetric water content %'

figure(2)
plot(m(:,1),m(:,2),'.')
hold on
xlabel 'temperature ^oC'
ylabel 'volumetric water content %'

T=m(:,1);%-5:0.1:1; %temperature in (C)
% 19cm 
a=0.06; b=-0.18; varnod=0.39; %volumetric soil moisture content
% 51cm a=0.06; b=-0.12; varnod=0.38; %volumetric soil moisture content
% 99cm  a=0.001; b=-0.01; varnod=0.4; %volumetric soil moisture content
t_f=-(varnod/a)^(1/b)
for i=1:length(T)
   if T(i)>t_f, theta(i)=varnod; %For temperatures greater than -0.01, UW function = vol. Moist. Cont.
   else theta(i)=a*abs(T(i))^b;  %For temperatures lower than -0.01, employ empirical law (Lovell 1957)
   end
end
%theta=theta/varnod; %Volumetric unfrozen water function
plot(T,theta,'-r','LineWidth',2); %plotting the graph
hold off
xlabel( 'temperature'); ylabel( 'volumetric UWC')
title(['a=',strcat(num2str(a)),', b=',strcat(num2str(b)),', VWC=',strcat(num2str(varnod))]);
grid on
%axis([-1.5 0.5 0 0.45])
