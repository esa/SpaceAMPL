%Programmed by:         Dario Izzo
%                       Advanced Concepts Team
%                       Advanced Concepts Team stagiaire
%Date:                  11/2007
%Revision:              1
%Tested by:             ----------
%
%
%
%
%Usage:     h=GTOC3_eph(ast,MJD2000)
%
%
%X-Ref :    needs the workspace contained in the file gtoc3.mat to
%know the asteroids elements



function [R,V]=gtoc3_eph(AstRow,MJD2000)

global AsteroidData
%JPL assigned constants
muSUN = 1.32712440018e11;
AU = 1.49597870691e8;


%This finds the ephemerides
epoch=AsteroidData(AstRow,8);
a=AsteroidData(AstRow,2)*AU;
e=AsteroidData(AstRow,3);
i=AsteroidData(AstRow,4);
W=AsteroidData(AstRow,5);
w=AsteroidData(AstRow,6);
M0=AsteroidData(AstRow,7);
epoch=epoch-51544; %CHECK IF CORRECT FROM JPL PROBLEM DESCRIPTION this should convert MJD to MJD2000
DT=(MJD2000-epoch)*60*60*24;
n=sqrt(muSUN/a^3);
M0=M0/180*pi;
M=M0+n*DT;
M=mod(M,2*pi);
E=M2E(M,e);
i = i/180*pi;
W = W/180*pi;
w = w/180*pi;
theta = 2 * atan( sqrt( (1+e)/(1-e)) * tan(E/2) );
gamma = atan ((e*sin(theta))/(1 + e*cos(theta)));
r = a * (1-e^2) / (1 + e *cos(theta));
v = sqrt (2*muSUN/r - muSUN / a);
x = r *(cos(theta + w)*cos(W)-sin(theta+w)*cos(i)*sin(W));
y = r *(cos(theta+w)*sin(W)+sin(theta+w)*cos(i)*cos(W));
z = r *(sin(theta+w)*sin(i));
vx = v * (-sin(theta+w-gamma)*cos(W)-cos(theta+w-gamma)*cos(i)*sin(W));
vy = v * (-sin(theta+w-gamma)*sin(W)+cos(theta+w-gamma)*cos(i)*cos(W));
vz = v * (cos(theta+w-gamma)*sin(i));

R = [x,y,z]';
V = [vx,vy,vz]';



