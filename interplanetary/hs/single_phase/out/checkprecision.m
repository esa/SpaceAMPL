%This script was written for the gtoc3 competition and loads the results
%from the ampl optimisation (Solution.out), the planet ephemerides
%(PlaParam.out), the units (units.out) and the trajectory times (Times.out),
%and integrate forward in time using Runge Kutta the equations of motion to
%check the numerical precision of the solution
%
%NOTE: tf is in MJD2000
function [TOUT,YOUT]=checkprecision(ast,tf)


ast1=ast(1);
ast2=ast(2);
Solution = load('Solution.out');
param = load('PlaParam.out');
Times = load('Times.out');
units = load('units.out');
load gtoc3.mat;

%set the units
R=1.49597870691E8; 		    % Length unit (in km)
MU=1.32712440018E11;		% Grav. parameter unit(km^3/sec^2)
V=sqrt(MU/R);			    % km/sec
A=V^2/R;				    % km/sec^2
M=units(3);
F=A*M;

%set the times
t0 = Times(1);
if nargin==1
    tf = Times(2) + t0;
    Tt  = Times(2) * 60*60*24; %simulation time
    Ttt=Tt;                    %times contained in the solution file
else
    Ttt  = Times(2) * 60*60*24; %simulation time
    Tt = (tf-t0) * 60*60*24;    %times contained in the solution file
end



%set nodes and collocation times
tcontrols = linspace(0,Ttt,size(Solution,1));

%set the controls
ux = Solution(:,8)*F;
uy = Solution(:,9)*F;
uz = Solution(:,10)*F;

%set initial conditions and TSPAN (the times at wich the solution is
%requested
y0=[Solution(1,1:3)*R Solution(1,4:6)*V Solution(1,7)*M];
TSPAN=[0:24*60*60:Tt];
if TSPAN(end)==Tt
else
    TSPAN = [TSPAN Tt];
end

%calls the integrator
options = odeset('RelTol',1e-9,'AbsTol',1e-9);
[TOUT,YOUT] = ode45(@(t,y) cartdynamic(t,y,tcontrols,ux,uy,uz),TSPAN,y0,options);

ur=sqrt((ux.^2+uy.^2+uz.^2));
cosx=ux./ur;
cosy=uy./ur;
cosz=uz./ur;
phi=angle(cosx+j*cosy);
uri = interp1(tcontrols,ur,TOUT);
phii = interp1(tcontrols, phi, TOUT);
thetai = interp1(tcontrols,acos(cosz),TOUT);
uxi = uri.*sin(thetai).*cos(phii);
uyi = uri.*sin(thetai).*sin(phii);
uzi = uri.*cos(thetai);

%Visualize the result
plot3(YOUT(:,1),YOUT(:,2),YOUT(:,3),'k')
hold
param(1,1)=param(1,1)*R;
param(2,1)=param(2,1)*R;
orbitpar(param(1,1:6));
orbitpar(param(2,1:6));
plot3(Solution(:,1)*R,Solution(:,2)*R,Solution(:,3)*R,'r')

%Evaluates the violation
if nargin == 1
    [R0,V0]=gtoc3_Eph(ast1,t0);
    [R1,V1]=gtoc3_Eph(ast2,tf);
    errR0=norm(R0'-YOUT(1,1:3))
    errV0=norm(V0'-YOUT(1,4:6))
    errRf=norm(R1'-YOUT(end,1:3))
    errVf=norm(V1'-YOUT(end,4:6))
else
    [R0,V0]=gtoc3_Eph(ast1,t0);
    errR0=norm(R0'-YOUT(1,1:3))
    errV0=norm(V0'-YOUT(1,4:6))
    check=[YOUT(end,1:3)/R YOUT(end,4:6)/V YOUT(end,7)/M]';
    tI=tf;
    tT=Times(2)+t0-tf;
    sprintf('%s%17.16f;\n%s%17.16f;\n%s%17.16f;\n%s%17.16f;\n%s%17.16f;\n%s%17.16f;\n%s%17.16f;\n%s%17.16f%s;\n%s%17.16f%s;'...
        ,'param x0 :=',check(1),'param y0 :=',check(2),'param z0 :=',check(3),'param dx0 :=', ...
        check(4),'param dy0 :=',check(5),'param dz0 :=',check(6),'param m0 :=',check(7),...
        'let timod :=',tf,' * d2u * f','let tfmod :=',Times(2)+t0-tf,' * d2u * f')
end
TOUT=TSPAN'/24/60/60+Times(1);
writesol=[TSPAN'/24/60/60+mjd20002mjd(Times(1)) YOUT uxi uyi uzi];
save SOLUTION.txt writesol -ascii -double -tabs


function dy=cartdynamic(t,y,tc,ux,uy,uz)

Isp=3000;
g0=9.80665*1e-3;
MU=1.32712440018E11;

ur=sqrt((ux.^2+uy.^2+uz.^2));
cosx=ux./ur;
cosy=uy./ur;
cosz=uz./ur;
phi=angle(cosx+j*cosy);
uri = interp1(tc,ur,t);
phii = interp1(tc, phi, t);
thetai = interp1(tc,acos(cosz),t);
uxi = uri.*sin(thetai).*cos(phii);
uyi = uri.*sin(thetai).*sin(phii);
uzi = uri.*cos(thetai);

dy(1) = y(4);
dy(2) = y(5);
dy(3) = y(6);
dy(4) = - y(1) * MU / (y(1)^2+y(2)^2+y(3)^2)^(3/2) + uxi/y(7);
dy(5) = - y(2) * MU / (y(1)^2+y(2)^2+y(3)^2)^(3/2) + uyi/y(7);
dy(6) = - y(3) * MU / (y(1)^2+y(2)^2+y(3)^2)^(3/2) + uzi/y(7);
dy(7) = - sqrt(uxi^2+uyi^2+uzi^2)/Isp/g0;
dy=dy';

