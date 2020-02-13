%{
some context: This was written for an assignment to find dynamic stability
variables for a thin panel with finite element modeling.
%}

clear; close all;

mTL=[156 -22;-22 4];
mTR=[54 13;-13 3];
mBL=[54 -13;13 -3];
mBR=[156 22;22 4];

mIN=mTL+mBR;

kTL=[6 -3;-3 2];
kTR=[-6 -3;3 1];
kBL=[-6 3;-3 1];
kBR=[6 3;3 2];

kIN=kTL+kBR;

ZERO=[0 0;0 0];

syms s

NT=[1-3*s^2+2*s^3;-s+2*s^2-s^3;3*s^2-2*s^3;s^2-s^3];
A=[-6*s+6*s^2, -1+4*s-3*s^2, 6*s-6*s^2, 2*s-3*s^2];
NTA=NT*A;

DefInt=int(NTA,s,0,1);

NTA_TL=[-1/2 -1/10;1/10 0];
NTA_TR=[1/2 1/10;-1/10 -1/60];
NTA_BL=[-1/2 1/10;-1/10 1/60];
NTA_BR=[1/2 -1/10;1/10 0];

NTA_IN=NTA_TL+NTA_BR;

%Three elements...
KB3=2*3^3*[kTL kTR ZERO ZERO;kBL kIN kTR ZERO;ZERO kBL kIN kTR;ZERO ZERO kBL kBR];
Ka3=[NTA_TL NTA_TR ZERO ZERO;NTA_BL NTA_IN NTA_TR ZERO;ZERO NTA_BL NTA_IN NTA_TR;ZERO ZERO NTA_BL NTA_BR];
M3=(1/(3*420))*[mTL mTR ZERO ZERO;mBL mIN mTR ZERO;ZERO mBL mIN mTR;ZERO ZERO mBL mBR];

%Boundary conditions applied...
KB3=[KB3(2:6,2:6) KB3(2:6,8);KB3(8,2:6) KB3(8,8)];
Ka3=[Ka3(2:6,2:6) Ka3(2:6,8);Ka3(8,2:6) Ka3(8,8)];
M3=[M3(2:6,2:6) M3(2:6,8);M3(8,2:6),M3(8,8)];

%Part d: Kq=qKaq
[V,D]=eig(KB3,Ka3);
imag(D)

N=600;

for i=0:N
    q(i+1)=i;
    Keff=KB3-q(i+1)*Ka3;
    [V,D]=eig(Keff,M3);
   
    sig5(i+1)=real(D(5,5));
    sig6(i+1)=real(D(6,6));
    
    w5(i+1)=imag(D(5,5));
    w6(i+1)=imag(D(6,6));
end


%At first, all of them were plotted to determine which were the smallest.

figure(1)
hold on
plot(q,sig5)
plot(q,sig6)
title('[3 element] \sigma vs. q')
xlabel('q')
ylabel('\sigma')
legend('\sigma_5','\sigma_6')

figure(2)
hold on
plot(q,w5)
plot(q,w6)
title('[3 element] \omega vs. q')
xlabel('q')
ylabel('\sigma')
legend('\omega_5','\omega_6')
