%{
some context: This code was written for a lab involving a hypersonic wind
tunnel.
%}

%{
XCQ-16:  10.035 mV/psi
XCQ-17:  10.075 mV/psi
XCE-121:  10.025 mV/psi
XCE-122:  9.949 mV/psi
%}

clear;close all;
load('C:\Users\redli\Documents\sko0l\ENAE464 - Aero Lab II\L6\ENAE464_PitotRake_Test2.mat')

%Sensitivity
sQ16=10.035;
sQ17=10.075;
sE121=10.025;
sE122=9.949;

p1=P1/6895;
p4=P4/6895;

%Find trigger
[M,ind]=max(trigger);
timestart=time(ind);
%timeend=time(ind+700);
%This was just a check
ind=ind+100;
%Cut to "relevant" data
timecut=time(ind:ind+700);

cut16=XCQ16(ind:ind+700);
cut17=XCQ17(ind:ind+700);
cut121=XCE121(ind:ind+700);
cut122=XCE122(ind:ind+700);

cutR1=Reservoir1(ind+49:ind+749);
cutR2=Reservoir2(ind+49:ind+749);
%accounting for the 0.49ms shift

%Conversion from V to psi
p16=1000*cut16/sQ16;
p17=1000*cut17/sQ17;
p121=1000*cut121/sE121;
p122=1000*cut122/sE122;

res1=1000*cutR1/10;
res2=1000*cutR2/10;

%Plots
figure(1)
hold on
plot(timecut,p16)
plot(timecut,p17)
plot(timecut,p121)
plot(timecut,p122)
legend('XCQ16','XCQ17','XCE121','XCE122')
title('Test Section')
ylabel('Pressure (psi)')
xlabel('Time (sec)')

figure(2)
hold on
plot(timecut,res1)
plot(timecut,res2)
legend('Res1','Res2')
title('Reservoir')
ylabel('Pressure (psi)')
xlabel('Time (sec)')

%Experimental Data from Pitot Tubes
avg16=mean(p16(100:end));
avg17=mean(p17(100:end));
avg121=mean(p121(100:end));
avg122=mean(p122(100:end));

avgR1=mean(res1(100:end-200));
avgR2=mean(res2(100:end-200));
avgRes=mean([avgR1,avgR2]);

xP02=mean([avg16,avg17,avg121,avg122]);

y=1.4;
y1=y;
M1=6;

%Rayleigh Pitot Formula
RP=((((y+1)*M1)^2/(4*y*M1^2-2*(y-1)))^(y/(y-1)))*((1-y+2*y*M1^2)/(y+1));
xP1=xP02/RP;

xP01=xP1*(1+(y-1)*M1^2/2)^(y/(y-1));

M2=sqrt((1+(y-1)*M1^2/2)/(y*M1^2-(y-1)/2));

xP2=xP1*xP02/(xP01*((1+(y-1)*M2^2/2)/(1+(y-1)*M1^2/2))^(y/(y-1)));
%xP01 "should" match P5


%Theoretical Calculations using P4/P1
Ratio4=p4/p1;
y4=gamma4;

%Speed of sound
R1=1717;
R4=12421;
RatioA=sqrt(y1*R1)/sqrt(y4*R4);
%Temperature cancels out

Func=@thing;
x=fsolve(Func,4);
%x = P2/P1
Ratio2=x;
tP2=Ratio2*p1;

Ms=sqrt((y+1)/(2*y)*(Ratio2-1)+1);

tP5=p1*((2*y1*Ms^2-(y1-1))/(y1+1))*(((3*y1-1)*Ms^2-2*(y1-1))/((y1-1)*Ms^2+2));

%Matching M1
M1g=7.5;
RPg=((((y+1)*M1g)^2/(4*y*M1g^2-2*(y-1)))^(y/(y-1)))*((1-y+2*y*M1g^2)/(y+1));
gP1=xP02/RPg;
gP01=gP1*(1+(y-1)*M1g^2/2)^(y/(y-1));

function F=thing(X)
y1=1.4;
y4=1.65;
RatioA=0.3425;
Ratio4=26.4814;
F=(X*(1-((y4-1)*RatioA*(X-1))/sqrt(2*y1*(2*y1+(y1+1)*(X-1))))^(-2*y4/(y4-1)))-Ratio4;
end





