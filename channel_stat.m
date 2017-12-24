clc
clear all
close all
A=load('/home/om1014/Project/turbulent/channel/upert.dat');
Y_plus=A(:,1);
U_stat=A(:,2);
UU_stat=A(:,3);
A=load('/home/om1014/Project/turbulent/channel/vpert.dat');
V_stat=A(:,2);

A=load('/home/om1014/Project/turbulent/channel/wpert.dat');
W_stat=A(:,2);


% A(abs(A)>=10000)=0;
plot(U_stat,Y_plus,'black','linewidth',2)

hold on
plot(V_stat,Y_plus,'r','linewidth',2)
plot(W_stat,Y_plus,'b','linewidth',2)
axis([0 max(U_stat)+0.5 0 80])
legend ('u^{,2}_{rms}','v^{,2}_{rms}','w^{,2}_{rms}')
[max(U_stat) max(V_stat) max(W_stat)]