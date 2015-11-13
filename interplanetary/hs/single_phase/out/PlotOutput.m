%This matlab script is associated to the AMPL model Pl2PlCartesian.mod and
%essentially plots the output solution. It needs SFToolbox V1.o or later
%(in particular it uses the orbit plotting routines)

initASTRO;
param=load('PlaParam.out');
units=load('units.out');

x=load('Solution.out');
times=load('Times.out');
xG=load('InitialGuess.out');
timesG=load('TimesGuess.out');

t=linspace(times(1),times(1)+times(2),size(x,1));
tG=linspace(timesG(1),timesG(1)+timesG(2),size(xG,1));

%read non dimensional units
R=units(1);
V=units(2);
M=units(3);
T=R/V;
A=V/T;
F=A*M;


%read ad plot Initial Guess
subplot(2,2,1);
plot3(xG(:,1),xG(:,2),xG(:,3),'r');
hold on
quiver3(xG(:,1),xG(:,2),xG(:,3),xG(:,8),xG(:,9),xG(:,10));

orbitpar(param(1,1:6));
orbitpar(param(2,1:6));
param(1,1)=param(1,1)*AU;
param(2,1)=param(2,1)*AU;
pospar(param(1,1:6),param(1,7),timesG(1),mu(11),1);
pospar(param(2,1:6),param(2,7),timesG(1)+timesG(2),mu(11),2);
axis square
title('Initial Guess Traj')

subplot(2,2,2);
param(1,1)=param(1,1)/AU;
param(2,1)=param(2,1)/AU;
plot3(xG(:,1),xG(:,2),xG(:,3),'r');
hold on
orbitpar(param(1,1:6));
orbitpar(param(2,1:6));
param(1,1)=param(1,1)*AU;
param(2,1)=param(2,1)*AU;
pospar(param(1,1:6),param(1,7),timesG(1),mu(11),1);
pospar(param(2,1:6),param(2,7),timesG(1)+timesG(2),mu(11),2);
axis square
title('Initial Guess Traj No Thrust')



subplot(2,2,3)
plot(tG,sqrt(xG(:,8).^2+xG(:,9).^2+xG(:,10).^2)*F*1000)
title('Initial Guess Thrust')
xlim([timesG(1) timesG(1)+timesG(2)])
subplot(2,2,4)
plot(tG,xG(:,7)*M)
title('Initial Guess Mass')
xlim([timesG(1) timesG(1)+timesG(2)])
param(1,1)=param(1,1)/AU;
param(2,1)=param(2,1)/AU;

figure
subplot(2,2,1)
%Plot Solution

plot3(x(:,1),x(:,2),x(:,3),'r');
hold on
quiver3(x(:,1),x(:,2),x(:,3),x(:,8),x(:,9),x(:,10));
orbitpar(param(1,1:6));
orbitpar(param(2,1:6));
param(1,1)=param(1,1)*AU;
param(2,1)=param(2,1)*AU;
pospar(param(1,1:6),param(1,7),times(1),mu(11),1);
pospar(param(2,1:6),param(2,7),times(1)+times(2),mu(11),2);
axis square
title('Solution')


subplot(2,2,2)
param(1,1)=param(1,1)/AU;
param(2,1)=param(2,1)/AU;
plot3(x(:,1),x(:,2),x(:,3),'r');
hold on
orbitpar(param(1,1:6));
orbitpar(param(2,1:6));
param(1,1)=param(1,1)*AU;
param(2,1)=param(2,1)*AU;
pospar(param(1,1:6),param(1,7),times(1),mu(11),1);
pospar(param(2,1:6),param(2,7),times(1)+times(2),mu(11),2);
axis square
title('Solution No Thrust')


subplot(2,2,3)
plot(t,sqrt(x(:,8).^2+x(:,9).^2+x(:,10).^2)*F*1000)
title('Thrust')
xlim([times(1) times(1)+times(2)])
subplot(2,2,4)
plot(t,x(:,7)*M)
title('Mass')
xlim([times(1) times(1)+times(2)])

