%{
some context: This was written for an assignment to compare different time
steps for numerical integration of an undamped free vibration problem
(trapezoidal rule).
%}

clear; close all;

M=1;
C=0;
K=4*pi^2;
F=0;

dt=0.01;

time1=0:dt:8;
N=8/dt+1;
q=zeros(1,N);
qd=zeros(1,N);
qdd=zeros(1,N);

q(1)=0;
qd(1)=1;
qdd(1)=(F-C*qd(1)-K*q(1))/M;

for n=1:N-1
    q(n+1)=(F+C*((2/dt)*q(n)+qd(n))+M*(((2/dt)^2)*q(n)+(4/dt)*qd(n)+qdd(n)))/(K+(2/dt)*C+((2/dt)^2)*M);
    qd(n+1)=(2/dt)*(q(n+1)-q(n))-qd(n);
    qdd(n+1)=((2/dt)^2)*(q(n+1)-q(n))-(4/dt)*qd(n)-qdd(n);
end

exact=(1/(2*pi))*sin(2*pi*time1);
exvel=cos(2*pi*time1); 
exactE=0.5*K*exact.^2+0.5*M*exvel.^2;
nde1=(0.5*K*q.^2+0.5*M*qd.^2-exactE)./exactE;
q1=q;

%

dt=0.05;

time2=0:dt:8;
N=8/dt+1;
q=zeros(1,N);
qd=zeros(1,N);
qdd=zeros(1,N);

q(1)=0;
qd(1)=1;
qdd(1)=(F-C*qd(1)-K*q(1))/M;

for n=1:N-1
    q(n+1)=(F+C*((2/dt)*q(n)+qd(n))+M*(((2/dt)^2)*q(n)+(4/dt)*qd(n)+qdd(n)))/(K+(2/dt)*C+((2/dt)^2)*M);
    qd(n+1)=(2/dt)*(q(n+1)-q(n))-qd(n);
    qdd(n+1)=((2/dt)^2)*(q(n+1)-q(n))-(4/dt)*qd(n)-qdd(n);
end

exact2=(1/(2*pi))*sin(2*pi*time2);
exvel2=cos(2*pi*time2); 
exactE2=0.5*K*exact2.^2+0.5*M*exvel2.^2;
nde2=(0.5*K*q.^2+0.5*M*qd.^2-exactE2)./exactE2;
q2=q;

%

dt=0.1;

time3=0:dt:8;
N=8/dt+1;
q=zeros(1,N);
qd=zeros(1,N);
qdd=zeros(1,N);

q(1)=0;
qd(1)=1;
qdd(1)=(F-C*qd(1)-K*q(1))/M;

for n=1:N-1
    q(n+1)=(F+C*((2/dt)*q(n)+qd(n))+M*(((2/dt)^2)*q(n)+(4/dt)*qd(n)+qdd(n)))/(K+(2/dt)*C+((2/dt)^2)*M);
    qd(n+1)=(2/dt)*(q(n+1)-q(n))-qd(n);
    qdd(n+1)=((2/dt)^2)*(q(n+1)-q(n))-(4/dt)*qd(n)-qdd(n);
end

exact3=(1/(2*pi))*sin(2*pi*time3);
exvel3=cos(2*pi*time3); 
exactE3=0.5*K*exact3.^2+0.5*M*exvel3.^2;
nde3=(0.5*K*q.^2+0.5*M*qd.^2-exactE3)./exactE3;
q3=q;

%

dt=0.5;

time4=0:dt:8;
N=8/dt+1;
q=zeros(1,N);
qd=zeros(1,N);
qdd=zeros(1,N);

q(1)=0;
qd(1)=1;
qdd(1)=(F-C*qd(1)-K*q(1))/M;


for n=1:N-1
    q(n+1)=(F+C*((2/dt)*q(n)+qd(n))+M*(((2/dt)^2)*q(n)+(4/dt)*qd(n)+qdd(n)))/(K+(2/dt)*C+((2/dt)^2)*M);
    qd(n+1)=(2/dt)*(q(n+1)-q(n))-qd(n);
    qdd(n+1)=((2/dt)^2)*(q(n+1)-q(n))-(4/dt)*qd(n)-qdd(n);
end

exact4=(1/(2*pi))*sin(2*pi*time4);
exvel4=cos(2*pi*time4); 
exactE4=0.5*K*exact4.^2+0.5*M*exvel4.^2;
nde4=(0.5*K*q.^2+0.5*M*qd.^2-exactE4)./exactE4;
q4=q;



figure(1)
hold on
plot(time1,exact)
plot(time1,q1)
plot(time2,q2)
plot(time3,q3)
plot(time4,q4)

xlabel('Time (sec)')
ylabel('Displacement')
title('Trapezoidal Rule')
legend('exact','dt=0.01','dt=0.05','dt=0.1','dt=0.5')

figure(2)
hold on
plot(time1,nde1)
plot(time2,nde2)
plot(time3,nde3)
plot(time4,nde4)

xlabel('Time (sec)')
ylabel('Total Energy')
title('Trapezoidal Rule - Energy Conservation')
legend('dt=0.01','dt=0.05','dt=0.1','dt=0.5')

    
    