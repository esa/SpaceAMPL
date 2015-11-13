let {i in 2..n-1} ux[i]:=Tmax*0.0000001;
let {i in 2..n-1} uy[i]:=Tmax*0.0000001;
let {i in 2..n-1} uz[i]:=Tmax*0.0000001;

#--------------------------------------------------------------------------
#Initial Guess for the DE variables
let  {i in J} DE[i] := DM1;
#-----------------------------------------------------------------------

let ux[2]:=dx2/v2*Tmax/2* tf/(n-1);
let uy[2]:=dy2/v2*Tmax/2* tf/(n-1);
let uz[2]:=dz2/v2*Tmax/2* tf/(n-1);
let ux[3]:=dx3/v3*Tmax/2* tf/(n-1);
let uy[3]:=dy3/v3*Tmax/2* tf/(n-1);
let uz[3]:=dz3/v3*Tmax/2* tf/(n-1);
let ux[4]:=dx4/v4*Tmax/2* tf/(n-1);
let uy[4]:=dy4/v4*Tmax/2* tf/(n-1);
let uz[4]:=dz4/v4*Tmax/2* tf/(n-1);
let ux[5]:=dx5/v5*Tmax/2* tf/(n-1);
let uy[5]:=dy5/v5*Tmax/2* tf/(n-1);
let uz[5]:=dz5/v5*Tmax/2* tf/(n-1);
let ux[6]:=dx6/v6*Tmax/2* tf/(n-1);
let uy[6]:=dy6/v6*Tmax/2* tf/(n-1);
let uz[6]:=dz6/v6*Tmax/2* tf/(n-1);
let ux[7]:=dx7/v7*Tmax/2* tf/(n-1);
let uy[7]:=dy7/v7*Tmax/2* tf/(n-1);
let uz[7]:=dz7/v7*Tmax/2* tf/(n-1);
let ux[8]:=dx8/v8*Tmax/2* tf/(n-1);
let uy[8]:=dy8/v8*Tmax/2* tf/(n-1);
let uz[8]:=dz8/v8*Tmax/2* tf/(n-1);
let ux[9]:=dx9/v9*Tmax/2* tf/(n-1);
let uy[9]:=dy9/v9*Tmax/2* tf/(n-1);
let uz[9]:=dz9/v9*Tmax/2* tf/(n-1);
let ux[10]:=dx10/v10*Tmax/2* tf/(n-1);
let uy[10]:=dy10/v10*Tmax/2* tf/(n-1);
let uz[10]:=dz10/v10*Tmax/2* tf/(n-1);
let ux[11]:=dx11/v11*Tmax/2* tf/(n-1);
let uy[11]:=dy11/v11*Tmax/2* tf/(n-1);
let uz[11]:=dz11/v11*Tmax/2* tf/(n-1);
let ux[12]:=dx12/v12*Tmax/2* tf/(n-1);
let uy[12]:=dy12/v12*Tmax/2* tf/(n-1);
let uz[12]:=dz12/v12*Tmax/2* tf/(n-1);
let ux[13]:=dx13/v13*Tmax/2* tf/(n-1);
let uy[13]:=dy13/v13*Tmax/2* tf/(n-1);
let uz[13]:=dz13/v13*Tmax/2* tf/(n-1);
let ux[14]:=dx14/v14*Tmax/2* tf/(n-1);
let uy[14]:=dy14/v14*Tmax/2* tf/(n-1);
let uz[14]:=dz14/v14*Tmax/2* tf/(n-1);
let ux[15]:=dx15/v15*Tmax/2* tf/(n-1);
let uy[15]:=dy15/v15*Tmax/2* tf/(n-1);
let uz[15]:=dz15/v15*Tmax/2* tf/(n-1);
let ux[16]:=dx16/v16*Tmax/2* tf/(n-1);
let uy[16]:=dy16/v16*Tmax/2* tf/(n-1);
let uz[16]:=dz16/v16*Tmax/2* tf/(n-1);
let ux[17]:=dx17/v17*Tmax/2* tf/(n-1);
let uy[17]:=dy17/v17*Tmax/2* tf/(n-1);
let uz[17]:=dz17/v17*Tmax/2* tf/(n-1);
let ux[18]:=dx18/v18*Tmax/2* tf/(n-1);
let uy[18]:=dy18/v18*Tmax/2* tf/(n-1);
let uz[18]:=dz18/v18*Tmax/2* tf/(n-1);
let ux[19]:=dx19/v19*Tmax/2* tf/(n-1);
let uy[19]:=dy19/v19*Tmax/2* tf/(n-1);
let uz[19]:=dz19/v19*Tmax/2* tf/(n-1);
let ux[20]:=dx20/v20*Tmax/2* tf/(n-1);
let uy[20]:=dy20/v20*Tmax/2* tf/(n-1);
let uz[20]:=dz20/v20*Tmax/2* tf/(n-1);
let ux[21]:=dx21/v21*Tmax/2* tf/(n-1);
let uy[21]:=dy21/v21*Tmax/2* tf/(n-1);
let uz[21]:=dz21/v21*Tmax/2* tf/(n-1);
let ux[22]:=dx22/v22*Tmax/2* tf/(n-1);
let uy[22]:=dy22/v22*Tmax/2* tf/(n-1);
let uz[22]:=dz22/v22*Tmax/2* tf/(n-1);
let ux[23]:=dx23/v23*Tmax/2* tf/(n-1);
let uy[23]:=dy23/v23*Tmax/2* tf/(n-1);
let uz[23]:=dz23/v23*Tmax/2* tf/(n-1);
let ux[24]:=dx24/v24*Tmax/2* tf/(n-1);
let uy[24]:=dy24/v24*Tmax/2* tf/(n-1);
let uz[24]:=dz24/v24*Tmax/2* tf/(n-1);
let ux[25]:=dx25/v25*Tmax/2* tf/(n-1);
let uy[25]:=dy25/v25*Tmax/2* tf/(n-1);
let uz[25]:=dz25/v25*Tmax/2* tf/(n-1);
let ux[26]:=dx26/v26*Tmax/2* tf/(n-1);
let uy[26]:=dy26/v26*Tmax/2* tf/(n-1);
let uz[26]:=dz26/v26*Tmax/2* tf/(n-1);
let ux[27]:=dx27/v27*Tmax/2* tf/(n-1);
let uy[27]:=dy27/v27*Tmax/2* tf/(n-1);
let uz[27]:=dz27/v27*Tmax/2* tf/(n-1);
let ux[28]:=dx28/v28*Tmax/2* tf/(n-1);
let uy[28]:=dy28/v28*Tmax/2* tf/(n-1);
let uz[28]:=dz28/v28*Tmax/2* tf/(n-1);
subject to
	thrustON{i in 2..n-1}: uT[i] <= Tmax*tf/(n-1);

minimize
	position: (xf-xn)^2+(yf-yn)^2+(zf-zn)^2+(dxf-dxn)^2+(dyf-dyn)^2+(dzf-dzn)^2;

drop FinalPositionx;
drop FinalPositiony;
drop FinalPositionz;
drop FinalVelocityx;
drop FinalVelocityy;
drop FinalVelocityz;
#--------------------------------------------------------------------------
solve;
#-----------------------------------------------------------------------
#--------------------------------------------------------------------------
#Print The Initial Guess x,y,z,dx,dy,dz,m,ux,uy,uz variables

param m{i in I} := 1;
include include/writeinitialguess.inc;
purge m;

#Print the initial and final times
printf "%17.16e, %17.16e \n", ti/d2u , tF-ti/d2u > out/TimesGuess.out;
close out/TimesGuess.out;
#------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Clean up
unfix timod;
unfix tfmod;
restore FinalPositionx;
restore FinalPositiony;
restore FinalPositionz;
restore FinalVelocityx;
restore FinalVelocityy;
restore FinalVelocityz;
drop thrustON;
drop position;
