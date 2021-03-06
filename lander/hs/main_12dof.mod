#####################################################################################
#
# Problem:  Optimal Control Problem (OCP)
# Dynamics: Spacecraft with four independent, symmetric nozzles
# Transcription: Hermite-Simpson
#
# Author: Dario Izzo (Mar 2021) 
#
#####################################################################################

#Sets---------------------------------------------------------
    set vI := 1..3;
    set vJ := 1..4;
#-------------------------------------------------------------

#Parameters---------------------------------------------------
#Generic
    param n default          50;                   # Number of nodes
    param g default          1.6229;               # [m/s^2] Gravitational Acceleration, default value is Moon
    param epsilon default    0.01;                 # Tunes the smoothing element
    param pi:=               4*atan(1);
    param Ispg0 default      3049;                 # [m/s] veff. product of the specific impulse with g0

#Spaceraft params   
    param Ixx default          30550 * 0.62;                 # [kg m^2] Moment of Inertia Tensor Component - xx
    param Iyy default          33014 * 0.62;                 # [kg m^2] Moment of Inertia Tensor Component - yy
    param Izz default          33826 * 0.62;                 # [kg m^2] Moment of Inertia Tensor Component - zz
    
    param L default            3.0;                          # [m] Distance thrust center of mass
    
    param maxthrust default           44000;                 # [N] Max thrust
    param max_secondary default       44;                    # [N] Max secondary thrust for the attitude control    
    
    param m0 default          9472.06;                       # [Kg] initial mass
    
#State constraint params
    param maxphi   default pi/6.0;
    param maxtheta default pi/6.0;
    param maxpsi   default pi/6.0;
    
#Initial conditions
    param x0 default           0.0;        # [m] Initial x
    param y0 default           0.0;        # [m] Initial y
    param z0 default           -2300;        # [m] Initial z
    
    param vx0 default          70.0;        # [m/s] Initial vx
    param vy0 default          0.0;        # [m/s] Initial vy
    param vz0 default          44.0;        # [m/s] Initial vz
    
    param phi0 default         0.0;        # [rad] Initial phi
    param theta0 default       0.0;        # [rad] Initial theta
    param psi0 default         0.0;        # [rad] Initial psi
    
    param p0 default           0.0;        # [rad] Initial phi rate
    param q0 default           0.0;        # [rad] Initial theta rate
    param r0 default           0.0;        # [rad] Initial psi rate
        
#Final conditions
    param xn default           0.0;        # [m] Final x
    param yn default           0.0;        # [m] Final y
    param zn default           0.0;        # [m] Final z
    
    param vxn default          0.0;        # [m/s] Final vx
    param vyn default          0.0;        # [m/s] Final vy
    param vzn default          0.0;        # [m/s] Final vz
    
    param phin default         0.0;       # [rad] Final phi
    param thetan default       0.0;        # [rad] Final theta
    param psin default         0.0;        # [rad] Final psi
    
    param pn default           0.0;        # [rad] Final phi rate
    param qn default           0.0;        # [rad] Final theta rate
    param rn default           0.0;        # [rad] Final psi rate
    
#Other
    param tn default          50.0;        #[s] Guess for the final time
#-------------------------------------------------------------

#Sets---------------------------------------------------------
    set I := 1..n;
    set J := 1..n-1;
#-------------------------------------------------------------

#Variables---------------------------------------------------
    var x {i in I};
    var y {i in I};
    var z {i in I}, <=0;
    
    var vx {i in I};
    var vy {i in I};
    var vz {i in I};
    
    var phi {i in I},   >= -maxphi,   <= maxphi;   # roll
    var theta {i in I}, >= -maxtheta, <= maxtheta; # pitch
    var psi {i in I},   >= -maxpsi, <= maxpsi;     # yaw
    
    var p {i in I};
    var q {i in I};
    var r {i in I};

    var m {i in I} := m0;
    
    var u  {i in I, j in vJ}, >=-1, <=1;  # 1 - main thruster, 2 - pitch, 3 - roll, 4 - yaw
    var um {i in J, j in vJ}, >=-1, <=1;  # 1 - main thruster, 2 - pitch, 3 - roll, 4 - yaw

#-------------------------------------------------------------

#Time variables-----------------------------------------------
    var tf, >=0;
    var dt = tf/(n-1);
    var timegrid{i in I} = dt*(i-1);
#-------------------------------------------------------------

#Objective----------------------------------------------------

        # For power, minimize Simpson's approximation to the integral:
        #
        #        \int{ f(t)dt }
        #     ~= \sum_{  dt/6 * f(t) + 4*f(t+dt/2)  + f(t+dt)  }
        #               for t=(dt,2*dt,3*dt...)
        #cost has the values at t = i*dt
        #costm has the values at t = i*dt + dt/2

    var cost{i in I}  = sum{j in vJ} u[i,j]^2;
    var costm{i in J} = sum{j in vJ} um[i,j]^2;
    var smoothing_term = dt/6 * sum{i in J} (cost[i]+4*costm[i]+cost[i+1]);

    minimize mass: -m[n] + epsilon*smoothing_term;
#-------------------------------------------------------------

#Dynamic at the grid points-----------------------------------
    var cp {i in I} = cos(phi[i]);
    var sp {i in I} = sin(phi[i]);
    var cq {i in I} = cos(theta[i]);
    var sq {i in I} = sin(theta[i]);
    var cr {i in I} = cos(psi[i]);
    var sr {i in I} = sin(psi[i]);
    
    var f1{i in I} = vx[i];
    var f2{i in I} = vy[i];
    var f3{i in I} = vz[i];
    
    var T{i in I}  = u[i, 1] * maxthrust / m[i];
    
    var f4{i in I} = -T[i]*(sp[i]*sr[i] + cp[i]*sq[i]*cr[i]);
    var f5{i in I} = -T[i]*(-sp[i]*cr[i] + cp[i]*sq[i]*sr[i]);
    var f6{i in I} = g-T[i]*(cp[i]*cq[i]);


    var f7{i in I} = p[i] + q[i]*(sp[i]*sq[i]/cq[i]) + r[i]*(cp[i]*sq[i]/cq[i]);
    var f8{i in I} = q[i]*cp[i] - r[i]*sp[i];
    var f9{i in I} = q[i]*sp[i]/cq[i] + r[i]*cp[i]/cq[i];
    
    var f10{i in I} = (q[i]*r[i]*(Iyy-Izz) + max_secondary * L * u[i,2]) / Ixx;
    var f11{i in I} = (p[i]*r[i]*(Izz-Ixx) + max_secondary * L * u[i,3]) / Iyy;
    var f12{i in I} = (p[i]*q[i]*(Ixx-Iyy) + max_secondary * L * u[i,4]) / Izz;

    var f13{i in I} = - (u[i, 1] * maxthrust + (u[i, 2] + u[i, 3] + u[i, 4]) * max_secondary) /  Ispg0;
#-----------------------------------------------------------------------

#State definition at mid-points via Hermite interpolation---------------
    var xm{i in J}      =   (    x[i] +     x[i+1])/2 + tf/(n-1)/8 * (f1[i] - f1[i+1]);
    var ym{i in J}      =   (    y[i] +     y[i+1])/2 + tf/(n-1)/8 * (f2[i] - f2[i+1]);
    var zm{i in J}      =   (    z[i] +     z[i+1])/2 + tf/(n-1)/8 * (f3[i] - f3[i+1]);
    
    var vxm{i in J}     =   (   vx[i] +    vx[i+1])/2 + tf/(n-1)/8 * (f4[i] - f4[i+1]);
    var vym{i in J}     =   (   vy[i] +    vy[i+1])/2 + tf/(n-1)/8 * (f5[i] - f5[i+1]);
    var vzm{i in J}     =   (   vz[i] +    vz[i+1])/2 + tf/(n-1)/8 * (f6[i] - f6[i+1]);
    
    var phim{i in J}    =   (phi[i]   + phi[i+1])/2   + tf/(n-1)/8 * (f7[i] - f7[i+1]);
    var thetam{i in J}  =   (theta[i] + theta[i+1])/2 + tf/(n-1)/8 * (f8[i] - f8[i+1]);
    var psim{i in J}    =   (psi[i]   + psi[i+1])/2   + tf/(n-1)/8 * (f9[i] - f9[i+1]);
    
    var pm{i in J}   =   (p[i] + p[i+1])/2   + tf/(n-1)/8 * (f10[i] - f10[i+1]);
    var qm{i in J}   =   (q[i] + q[i+1])/2   + tf/(n-1)/8 * (f11[i] - f11[i+1]);
    var rm{i in J}   =   (r[i] + r[i+1])/2   + tf/(n-1)/8 * (f12[i] - f12[i+1]);

    var mm{i in J}   =   (m[i] + m[i+1])/2   + tf/(n-1)/8 * (f13[i] - f13[i+1]);
    #-----------------------------------------------------------------------

#Dynamic at the mid-points----------------------------------------------
    var cpm {i in J} = cos(phim[i]);
    var spm {i in J} = sin(phim[i]);
    var cqm {i in J} = cos(thetam[i]);
    var sqm {i in J} = sin(thetam[i]);
    var crm {i in J} = cos(psim[i]);
    var srm {i in J} = sin(psim[i]);
    
    var f1m{i in J} = vxm[i];
    var f2m{i in J} = vym[i];
    var f3m{i in J} = vzm[i];
    
    var Tm{i in J}  =  um[i, 1] * maxthrust / mm[i];
    
    var f4m{i in J} = -Tm[i]*(spm[i]*srm[i] + cpm[i]*sqm[i]*crm[i]);   
    var f5m{i in J} = -Tm[i]*(-spm[i]*crm[i] + cpm[i]*sqm[i]*srm[i]);
    var f6m{i in J} = g-Tm[i]*(cpm[i]*cqm[i]);


    var f7m{i in J} = pm[i] + qm[i]*(spm[i]*sqm[i]/cqm[i]) + rm[i]*(cpm[i]*sqm[i]/cqm[i]);
    var f8m{i in J} = qm[i]*cpm[i] - rm[i]*spm[i];
    var f9m{i in J} = qm[i]*spm[i]/cqm[i] + rm[i]*cpm[i]/cqm[i];
    
    var f10m{i in J} = (qm[i]*rm[i]*(Iyy-Izz) + max_secondary * L * um[i,2]) / Ixx;
    var f11m{i in J} = (pm[i]*rm[i]*(Izz-Ixx) + max_secondary * L * um[i,3]) / Iyy;
    var f12m{i in J} = (pm[i]*qm[i]*(Ixx-Iyy) + max_secondary * L * um[i,4]) / Izz;

    var f13m{i in J} = - (um[i, 1] * maxthrust + (um[i, 2] + um[i, 3] + um[i, 4]) * max_secondary) /  Ispg0;
#-----------------------------------------------------------------------

#Simpson Formula---------------------------------------------------------
subject to
    dynamicx{i in J}:         x[i]  =    x[i+1] - tf/(n-1)/6*(f1[i] + f1[i+1] + 4*f1m[i]);
    dynamicy{i in J}:         y[i]  =    y[i+1] - tf/(n-1)/6*(f2[i] + f2[i+1] + 4*f2m[i]);
    dynamicz{i in J}:         z[i]  =    z[i+1] - tf/(n-1)/6*(f3[i] + f3[i+1] + 4*f3m[i]);
    
    dynamicvx{i in J}:       vx[i]  =    vx[i+1] - tf/(n-1)/6*(f4[i] + f4[i+1] + 4*f4m[i]);
    dynamicvy{i in J}:       vy[i]  =    vy[i+1] - tf/(n-1)/6*(f5[i] + f5[i+1] + 4*f5m[i]);
    dynamicvz{i in J}:       vz[i]  =    vz[i+1] - tf/(n-1)/6*(f6[i] + f6[i+1] + 4*f6m[i]);
    
    dynamicp{i in J}: phi[i]    = phi[i+1]   - tf/(n-1)/6*(f7[i] + f7[i+1] + 4*f7m[i]);
    dynamicq{i in J}: theta[i]  = theta[i+1] - tf/(n-1)/6*(f8[i] + f8[i+1] + 4*f8m[i]);
    dynamicr{i in J}: psi[i]    = psi[i+1]   - tf/(n-1)/6*(f9[i] + f9[i+1] + 4*f9m[i]);
    
    dynamicdp{i in J}: p[i]    = p[i+1]   - tf/(n-1)/6*(f10[i] + f10[i+1] + 4*f10m[i]);
    dynamicdq{i in J}: q[i]    = q[i+1]   - tf/(n-1)/6*(f11[i] + f11[i+1] + 4*f11m[i]);
    dynamicdr{i in J}: r[i]    = r[i+1]   - tf/(n-1)/6*(f12[i] + f12[i+1] + 4*f12m[i]);

    dynamicdm{i in J}: m[i]    = m[i+1]   - tf/(n-1)/6*(f13[i] + f13[i+1] + 4*f13m[i]);
#--------------------------------------------------------------------------

#Constraints------------------------------------------
    #Boundary Conditions
    
    #Initial
    subject to InitialPositionx :  x[1] = x0;
    subject to InitialPositiony :  y[1] = y0;
    subject to InitialPositionz :  z[1] = z0;
    
    subject to InitialVelocityx : vx[1] = vx0;
    subject to InitialVelocityy : vy[1] = vy0;
    subject to InitialVelocityz : vz[1] = vz0;
    
    subject to InitialPitch     :  phi[1]   = phi0;
    subject to InitialRoll      :  theta[1] = theta0;
    subject to InitialYaw       :  psi[1]   = psi0;
    
    subject to InitialPitchRate : p[1]   = p0;
    subject to InitialRollRate  : q[1]   = q0;
    subject to InitialYawRate   : r[1]   = r0;

    subject to InitialMass   : m[1]   = m0;

    
    #Final
    #subject to FinalPositionx :  x[n] = xn;
    #subject to FinalPositiony :  y[n] = yn;
    subject to FinalPositionz :  z[n] = zn;
    
    subject to FinalVelocityx : vx[n] = vxn;
    subject to FinalVelocityy : vy[n] = vyn;
    subject to FinalVelocityz : vz[n] = vzn;
    
    subject to FinalPitch     :  phi[n]   = phin;
    subject to FinalRoll      :  theta[n] = thetan;
    subject to FinalYaw       :  psi[n]   = psin;
    
    subject to FinalPitchRate : p[n]   = pn;
    subject to FinalRollRate  : q[n]   = qn;
    subject to FinalYawRate   : r[n]   = rn;

    #  main thruster only goes up
    subject to mainthruster{i in I}: u[i,1] >=0; # main thruster only goes up
    subject to mainthrusterm{i in J}: um[i,1] >=0; # main thruster only goes up
    
#-------------------------------------------------------------

#Guess-------------------------------------------------------
    let tf := tn;
    let {i in I, j in vJ}  u[i, j] := 0.5;
    let {i in J, j in vJ} um[i, j] := 0.5;
    let {i in I} m[i] := m0;

#-------------------------------------------------------------

#Solver Options-----------------------------------------------
    option solver snopt;
    option substout 0;
    option show_stats 1;
    options snopt_options "outlev=2 Major_iterations=1500 Superbasics=500";
#-------------------------------------------------------------
