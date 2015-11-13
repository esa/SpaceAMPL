%This matlab script is associated to the AMPL model FlyByPl2Pl and
%essentially plots the output solution. It needs SFToolbox V1.o or later
%(in particular it uses the orbit plotting routines). It reads from the
%current directory the files produced by the AMPL model output
%PlaParam.out, units.out, InitialGuess1.out TimesGuess1.out,
%InitialGuess2.out, TimesGuess2.out

initASTRO;
param=csvread('PlaParam.out');
units=csvread('units.out');

%read non dimensional units
R=units(1);
V=units(2);
M=units(3);
T=R/V;
A=V/T;
F=A*M;


%read ad plot Initial Guess 
xG=csvread('InitialGuess1.out');
timesG=csvread('TimesGuess1.out');

xG2=csvread('InitialGuess2.out');
timesG2=csvread('TimesGuess2.out');



subplot(2,2,1);
plot3(xG(:,1),xG(:,2),xG(:,3),'r');
hold on
quiver3(xG(:,1),xG(:,2),xG(:,3),xG(:,8),xG(:,9),xG(:,10));
plot3(xG2(:,1),xG2(:,2),xG2(:,3),'k');
quiver3(xG2(:,1),xG2(:,2),xG2(:,3),xG2(:,8),xG2(:,9),xG(:,10));


orbitpar(param(1,1:6));
orbitpar(param(2,1:6));
orbitpar(param(3,1:6));
param(1,1)=param(1,1)*AU;
param(2,1)=param(2,1)*AU;
param(3,1)=param(3,1)*AU;

pospar(param(1,1:6),param(1,7),timesG(1),mu(11),1);
pospar(param(2,1:6),param(2,7),timesG(1)+timesG(2),mu(11),2);

pospar(param(2,1:6),param(2,7),timesG2(1),mu(11),1);
pospar(param(3,1:6),param(3,7),timesG2(1)+timesG2(2),mu(11),2);

axis square
title('Initial Guess Traj')

subplot(2,2,2);
param(1,1)=param(1,1)/AU;
param(2,1)=param(2,1)/AU;
param(3,1)=param(3,1)/AU;

plot3(xG(:,1),xG(:,2),xG(:,3),'r');
hold on
quiver3(xG(:,1),xG(:,2),xG(:,3),xG(:,8),xG(:,9),xG(:,10));
plot3(xG2(:,1),xG2(:,2),xG2(:,3),'k');
quiver3(xG2(:,1),xG2(:,2),xG2(:,3),xG2(:,8),xG2(:,9),xG(:,10));


orbitpar(param(1,1:6));
orbitpar(param(2,1:6));
orbitpar(param(3,1:6));
param(1,1)=param(1,1)*AU;
param(2,1)=param(2,1)*AU;
param(3,1)=param(3,1)*AU;

pospar(param(1,1:6),param(1,7),timesG(1),mu(11),1);
pospar(param(2,1:6),param(2,7),timesG(1)+timesG(2),mu(11),2);

pospar(param(2,1:6),param(2,7),timesG2(1),mu(11),1);
pospar(param(3,1:6),param(3,7),timesG2(1)+timesG2(2),mu(11),2);
axis square
title('Initial Guess Traj XY')
view(0,90);


subplot(2,2,3)
plot(sqrt(xG(:,8).^2+xG(:,9).^2+xG(:,10).^2)*F*1000)
title('Initial Guess Thrust')
subplot(2,2,4)
plot(xG(:,7)*M)
title('Initial Guess Mass')

figure
subplot(2,2,2)
%Plot Solution
x1=csvread('Solution1.out');
param=csvread('PlaParam.out');
times1=csvread('Times1.out');

x2=csvread('Solution2.out');
times2=csvread('Times2.out');

plot3(x1(:,1),x1(:,2),x1(:,3),'r');
hold on
%quiver3(x1(:,1),x1(:,2),x1(:,3),x1(:,8),x1(:,9),x1(:,10));
plot3(x2(:,1),x2(:,2),x2(:,3),'k');
%quiver3(x2(:,1),x2(:,2),x2(:,3),x2(:,8),x2(:,9),x2(:,10));



orbitpar(param(1,1:6));
orbitpar(param(2,1:6));
orbitpar(param(3,1:6));
param(1,1)=param(1,1)*AU;
param(2,1)=param(2,1)*AU;
param(3,1)=param(3,1)*AU;
pospar(param(1,1:6),param(1,7),times1(1),mu(11),1);
pospar(param(2,1:6),param(2,7),times1(1)+times1(2),mu(11),2);

pospar(param(2,1:6),param(2,7),times2(1),mu(11),1);
pospar(param(3,1:6),param(3,7),times2(1)+times2(2),mu(11),2);
axis square
title('Solution')
 
subplot(2,2,3)
plot([sqrt(x1(:,8).^2+x1(:,9).^2+x1(:,10).^2); sqrt(x2(:,8).^2+x2(:,9).^2+x2(:,10).^2)] *F*1000)
title('Thrust')
subplot(2,2,4)
plot([x1(:,7); x2(:,7)]*M)
title('Mass')
 

