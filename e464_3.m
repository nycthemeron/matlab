%{
some context: This was written for a lab involving a pendulum on a rotary
beam. In this lab, we found parameters of the system and then used output
feedback control to stabilize the inverted pendulum.
%}

clear
close all
load('C:\Users\redli\Documents\ENAE464 - Aero Lab II\all_data.mat')

%%These three sets of data were tested one at a time.
%{
%%Pendulum only data:
%Trim the data so 4-5 periods are shown.
s1=10000;
e1=13000;

theta1=data1.theta(s1:e1);
time1=data1.time(s1:e1);
time1=time1-time1(1);

%x= [amplitude, frequency, phase, offset, decay-rate]
%These values are changed until the two plots line up. 
x01= [.16 3 3 0 .07];
out1=datafit(theta1,time1,x01);

omega1=out1(2)
sigma1=out1(5)

%%Rotary beam only data:
%Trim the data so 4-5 periods are shown.
s2=1;
e2=2500;
phi2=data2.phi(s2:e2);
time2=data2.time(s2:e2);
time2=time2-time2(1);


%These values are changed until the two plots line up.
x02= [.7 3 3 0 .2];
out2=datafit(phi2,time2,x02);

omega2=out2(2)
sigma2=out2(5)
%}
%%Rotary beam with disk weights data:
%Trim the data so 4-5 periods are shown.
s3=1;
e3=2500;
phi3=data3.phi(s3:e3);
time3=data3.time(s3:e3);
time3=time3-time3(1);

%These values are changed until the two plots line up.
x03= [.7 3 3 0 .2];
out3=datafit(phi3,time3,x03);

omega3=out3(2)
sigma3=out3(5)

function err = datafitfun(x,time,data,handle)
fitdata=x(1)*exp(-x(5)*time).*sin(x(2)*time+x(3))+x(4);
err=norm(fitdata-data);
end

function out = datafit(exp_data,time,x0)
plot(time,exp_data,'k.-','linewidth',2)
hold on, grid on
xlabel('time'), ylabel('data')
h=plot(time,zeros(length(exp_data)),'r','linewidth',2);
options=optimset('TolX',0.001);

out=fminsearch(@datafitfun,x0,options,time,exp_data,h);

fitdata = out(1)*exp(-out(5)*time).*sin(out(2)*time+out(3))+out(4);
set(h, 'ydata', fitdata)
drawnow
end