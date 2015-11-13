#--------------------------------------------------------------------------
# Node 1: Departure Node
var x1 =  x0;
var y1 =  y0;
var z1 =  z0;
var dx1 =  dx0 + ux[1];
var dy1 =  dy0 + uy[1];
var dz1 =  dz0 + uz[1];

#Basic definitions
var r1 = sqrt(x1^2+y1^2+z1^2);
var v1 = sqrt(dx1^2+dy1^2+dz1^2);
var a1 = 1 / (2/r1 - v1^2);
var sigma1 = x1*dx1+y1*dy1+z1*dz1;
var meanmotion1 = sqrt(1/a1^3);
var DM1 = meanmotion1 * dt;

#Lagrange Coefficients
var rvar1 = a1 + (r1-a1)*cos(DE[1]) + sigma1*sqrt(a1)*sin(DE[1]);
var F1 = 1 - a1/r1 * (1-cos(DE[1]));
var G1 = a1*sigma1*(1-cos(DE[1])) + r1*sqrt(a1)*sin(DE[1]);
var Ft1 = -sqrt(a1)/(r1*rvar1) * sin(DE[1]);
var Gt1 = 1 - a1/rvar1*(1-cos(DE[1]));

subject to KeplerEquations1: 
	DM1 - DE[1] - sigma1/sqrt(a1) * (1 - cos(DE[1])) + (1 - r1/a1)*sin(DE[1]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 2
var x2 = F1*x1 + G1*(dx1);
var y2 = F1*y1 + G1*(dy1);
var z2 = F1*z1 + G1*(dz1);
var dx2 = Ft1*x1 + Gt1*(dx1)+ux[2];
var dy2 = Ft1*y1 + Gt1*(dy1)+uy[2];
var dz2 = Ft1*z1 + Gt1*(dz1)+uz[2];

#Basic definitions
var r2 = sqrt(x2^2+y2^2+z2^2);
var v2 = sqrt(dx2^2+dy2^2+dz2^2);
var a2 = 1 / (2/r2 - v2^2);
var sigma2 = x2*dx2+y2*dy2+z2*dz2;
var meanmotion2 = sqrt(1/a2^3);
var DM2 = meanmotion2 * dt;

#Lagrange Coefficients
var rvar2 = a2 + (r2-a2)*cos(DE[2]) + sigma2*sqrt(a2)*sin(DE[2]);
var F2 = 1 - a2/r2 * (1-cos(DE[2]));
var G2 = a2*sigma2*(1-cos(DE[2])) + r2*sqrt(a2)*sin(DE[2]);
var Ft2 = -sqrt(a2)/(r2*rvar2) * sin(DE[2]);
var Gt2 = 1 - a2/rvar2*(1-cos(DE[2]));

subject to KeplerEquations2: 
	DM2 - DE[2] - sigma2/sqrt(a2) * (1 - cos(DE[2])) + (1 - r2/a2)*sin(DE[2]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 3
var x3 = F2*x2 + G2*(dx2);
var y3 = F2*y2 + G2*(dy2);
var z3 = F2*z2 + G2*(dz2);
var dx3 = Ft2*x2 + Gt2*(dx2)+ux[3];
var dy3 = Ft2*y2 + Gt2*(dy2)+uy[3];
var dz3 = Ft2*z2 + Gt2*(dz2)+uz[3];

#Basic definitions
var r3 = sqrt(x3^2+y3^2+z3^2);
var v3 = sqrt(dx3^2+dy3^2+dz3^2);
var a3 = 1 / (2/r3 - v3^2);
var sigma3 = x3*dx3+y3*dy3+z3*dz3;
var meanmotion3 = sqrt(1/a3^3);
var DM3 = meanmotion3 * dt;

#Lagrange Coefficients
var rvar3 = a3 + (r3-a3)*cos(DE[3]) + sigma3*sqrt(a3)*sin(DE[3]);
var F3 = 1 - a3/r3 * (1-cos(DE[3]));
var G3 = a3*sigma3*(1-cos(DE[3])) + r3*sqrt(a3)*sin(DE[3]);
var Ft3 = -sqrt(a3)/(r3*rvar3) * sin(DE[3]);
var Gt3 = 1 - a3/rvar3*(1-cos(DE[3]));

subject to KeplerEquations3: 
	DM3 - DE[3] - sigma3/sqrt(a3) * (1 - cos(DE[3])) + (1 - r3/a3)*sin(DE[3]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 4
var x4 = F3*x3 + G3*(dx3);
var y4 = F3*y3 + G3*(dy3);
var z4 = F3*z3 + G3*(dz3);
var dx4 = Ft3*x3 + Gt3*(dx3)+ux[4];
var dy4 = Ft3*y3 + Gt3*(dy3)+uy[4];
var dz4 = Ft3*z3 + Gt3*(dz3)+uz[4];

#Basic definitions
var r4 = sqrt(x4^2+y4^2+z4^2);
var v4 = sqrt(dx4^2+dy4^2+dz4^2);
var a4 = 1 / (2/r4 - v4^2);
var sigma4 = x4*dx4+y4*dy4+z4*dz4;
var meanmotion4 = sqrt(1/a4^3);
var DM4 = meanmotion4 * dt;

#Lagrange Coefficients
var rvar4 = a4 + (r4-a4)*cos(DE[4]) + sigma4*sqrt(a4)*sin(DE[4]);
var F4 = 1 - a4/r4 * (1-cos(DE[4]));
var G4 = a4*sigma4*(1-cos(DE[4])) + r4*sqrt(a4)*sin(DE[4]);
var Ft4 = -sqrt(a4)/(r4*rvar4) * sin(DE[4]);
var Gt4 = 1 - a4/rvar4*(1-cos(DE[4]));

subject to KeplerEquations4: 
	DM4 - DE[4] - sigma4/sqrt(a4) * (1 - cos(DE[4])) + (1 - r4/a4)*sin(DE[4]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 5
var x5 = F4*x4 + G4*(dx4);
var y5 = F4*y4 + G4*(dy4);
var z5 = F4*z4 + G4*(dz4);
var dx5 = Ft4*x4 + Gt4*(dx4)+ux[5];
var dy5 = Ft4*y4 + Gt4*(dy4)+uy[5];
var dz5 = Ft4*z4 + Gt4*(dz4)+uz[5];

#Basic definitions
var r5 = sqrt(x5^2+y5^2+z5^2);
var v5 = sqrt(dx5^2+dy5^2+dz5^2);
var a5 = 1 / (2/r5 - v5^2);
var sigma5 = x5*dx5+y5*dy5+z5*dz5;
var meanmotion5 = sqrt(1/a5^3);
var DM5 = meanmotion5 * dt;

#Lagrange Coefficients
var rvar5 = a5 + (r5-a5)*cos(DE[5]) + sigma5*sqrt(a5)*sin(DE[5]);
var F5 = 1 - a5/r5 * (1-cos(DE[5]));
var G5 = a5*sigma5*(1-cos(DE[5])) + r5*sqrt(a5)*sin(DE[5]);
var Ft5 = -sqrt(a5)/(r5*rvar5) * sin(DE[5]);
var Gt5 = 1 - a5/rvar5*(1-cos(DE[5]));

subject to KeplerEquations5: 
	DM5 - DE[5] - sigma5/sqrt(a5) * (1 - cos(DE[5])) + (1 - r5/a5)*sin(DE[5]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 6
var x6 = F5*x5 + G5*(dx5);
var y6 = F5*y5 + G5*(dy5);
var z6 = F5*z5 + G5*(dz5);
var dx6 = Ft5*x5 + Gt5*(dx5)+ux[6];
var dy6 = Ft5*y5 + Gt5*(dy5)+uy[6];
var dz6 = Ft5*z5 + Gt5*(dz5)+uz[6];

#Basic definitions
var r6 = sqrt(x6^2+y6^2+z6^2);
var v6 = sqrt(dx6^2+dy6^2+dz6^2);
var a6 = 1 / (2/r6 - v6^2);
var sigma6 = x6*dx6+y6*dy6+z6*dz6;
var meanmotion6 = sqrt(1/a6^3);
var DM6 = meanmotion6 * dt;

#Lagrange Coefficients
var rvar6 = a6 + (r6-a6)*cos(DE[6]) + sigma6*sqrt(a6)*sin(DE[6]);
var F6 = 1 - a6/r6 * (1-cos(DE[6]));
var G6 = a6*sigma6*(1-cos(DE[6])) + r6*sqrt(a6)*sin(DE[6]);
var Ft6 = -sqrt(a6)/(r6*rvar6) * sin(DE[6]);
var Gt6 = 1 - a6/rvar6*(1-cos(DE[6]));

subject to KeplerEquations6: 
	DM6 - DE[6] - sigma6/sqrt(a6) * (1 - cos(DE[6])) + (1 - r6/a6)*sin(DE[6]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 7
var x7 = F6*x6 + G6*(dx6);
var y7 = F6*y6 + G6*(dy6);
var z7 = F6*z6 + G6*(dz6);
var dx7 = Ft6*x6 + Gt6*(dx6)+ux[7];
var dy7 = Ft6*y6 + Gt6*(dy6)+uy[7];
var dz7 = Ft6*z6 + Gt6*(dz6)+uz[7];

#Basic definitions
var r7 = sqrt(x7^2+y7^2+z7^2);
var v7 = sqrt(dx7^2+dy7^2+dz7^2);
var a7 = 1 / (2/r7 - v7^2);
var sigma7 = x7*dx7+y7*dy7+z7*dz7;
var meanmotion7 = sqrt(1/a7^3);
var DM7 = meanmotion7 * dt;

#Lagrange Coefficients
var rvar7 = a7 + (r7-a7)*cos(DE[7]) + sigma7*sqrt(a7)*sin(DE[7]);
var F7 = 1 - a7/r7 * (1-cos(DE[7]));
var G7 = a7*sigma7*(1-cos(DE[7])) + r7*sqrt(a7)*sin(DE[7]);
var Ft7 = -sqrt(a7)/(r7*rvar7) * sin(DE[7]);
var Gt7 = 1 - a7/rvar7*(1-cos(DE[7]));

subject to KeplerEquations7: 
	DM7 - DE[7] - sigma7/sqrt(a7) * (1 - cos(DE[7])) + (1 - r7/a7)*sin(DE[7]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 8
var x8 = F7*x7 + G7*(dx7);
var y8 = F7*y7 + G7*(dy7);
var z8 = F7*z7 + G7*(dz7);
var dx8 = Ft7*x7 + Gt7*(dx7)+ux[8];
var dy8 = Ft7*y7 + Gt7*(dy7)+uy[8];
var dz8 = Ft7*z7 + Gt7*(dz7)+uz[8];

#Basic definitions
var r8 = sqrt(x8^2+y8^2+z8^2);
var v8 = sqrt(dx8^2+dy8^2+dz8^2);
var a8 = 1 / (2/r8 - v8^2);
var sigma8 = x8*dx8+y8*dy8+z8*dz8;
var meanmotion8 = sqrt(1/a8^3);
var DM8 = meanmotion8 * dt;

#Lagrange Coefficients
var rvar8 = a8 + (r8-a8)*cos(DE[8]) + sigma8*sqrt(a8)*sin(DE[8]);
var F8 = 1 - a8/r8 * (1-cos(DE[8]));
var G8 = a8*sigma8*(1-cos(DE[8])) + r8*sqrt(a8)*sin(DE[8]);
var Ft8 = -sqrt(a8)/(r8*rvar8) * sin(DE[8]);
var Gt8 = 1 - a8/rvar8*(1-cos(DE[8]));

subject to KeplerEquations8: 
	DM8 - DE[8] - sigma8/sqrt(a8) * (1 - cos(DE[8])) + (1 - r8/a8)*sin(DE[8]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 9
var x9 = F8*x8 + G8*(dx8);
var y9 = F8*y8 + G8*(dy8);
var z9 = F8*z8 + G8*(dz8);
var dx9 = Ft8*x8 + Gt8*(dx8)+ux[9];
var dy9 = Ft8*y8 + Gt8*(dy8)+uy[9];
var dz9 = Ft8*z8 + Gt8*(dz8)+uz[9];

#Basic definitions
var r9 = sqrt(x9^2+y9^2+z9^2);
var v9 = sqrt(dx9^2+dy9^2+dz9^2);
var a9 = 1 / (2/r9 - v9^2);
var sigma9 = x9*dx9+y9*dy9+z9*dz9;
var meanmotion9 = sqrt(1/a9^3);
var DM9 = meanmotion9 * dt;

#Lagrange Coefficients
var rvar9 = a9 + (r9-a9)*cos(DE[9]) + sigma9*sqrt(a9)*sin(DE[9]);
var F9 = 1 - a9/r9 * (1-cos(DE[9]));
var G9 = a9*sigma9*(1-cos(DE[9])) + r9*sqrt(a9)*sin(DE[9]);
var Ft9 = -sqrt(a9)/(r9*rvar9) * sin(DE[9]);
var Gt9 = 1 - a9/rvar9*(1-cos(DE[9]));

subject to KeplerEquations9: 
	DM9 - DE[9] - sigma9/sqrt(a9) * (1 - cos(DE[9])) + (1 - r9/a9)*sin(DE[9]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 10
var x10 = F9*x9 + G9*(dx9);
var y10 = F9*y9 + G9*(dy9);
var z10 = F9*z9 + G9*(dz9);
var dx10 = Ft9*x9 + Gt9*(dx9)+ux[10];
var dy10 = Ft9*y9 + Gt9*(dy9)+uy[10];
var dz10 = Ft9*z9 + Gt9*(dz9)+uz[10];

#Basic definitions
var r10 = sqrt(x10^2+y10^2+z10^2);
var v10 = sqrt(dx10^2+dy10^2+dz10^2);
var a10 = 1 / (2/r10 - v10^2);
var sigma10 = x10*dx10+y10*dy10+z10*dz10;
var meanmotion10 = sqrt(1/a10^3);
var DM10 = meanmotion10 * dt;

#Lagrange Coefficients
var rvar10 = a10 + (r10-a10)*cos(DE[10]) + sigma10*sqrt(a10)*sin(DE[10]);
var F10 = 1 - a10/r10 * (1-cos(DE[10]));
var G10 = a10*sigma10*(1-cos(DE[10])) + r10*sqrt(a10)*sin(DE[10]);
var Ft10 = -sqrt(a10)/(r10*rvar10) * sin(DE[10]);
var Gt10 = 1 - a10/rvar10*(1-cos(DE[10]));

subject to KeplerEquations10: 
	DM10 - DE[10] - sigma10/sqrt(a10) * (1 - cos(DE[10])) + (1 - r10/a10)*sin(DE[10]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 11
var x11 = F10*x10 + G10*(dx10);
var y11 = F10*y10 + G10*(dy10);
var z11 = F10*z10 + G10*(dz10);
var dx11 = Ft10*x10 + Gt10*(dx10)+ux[11];
var dy11 = Ft10*y10 + Gt10*(dy10)+uy[11];
var dz11 = Ft10*z10 + Gt10*(dz10)+uz[11];

#Basic definitions
var r11 = sqrt(x11^2+y11^2+z11^2);
var v11 = sqrt(dx11^2+dy11^2+dz11^2);
var a11 = 1 / (2/r11 - v11^2);
var sigma11 = x11*dx11+y11*dy11+z11*dz11;
var meanmotion11 = sqrt(1/a11^3);
var DM11 = meanmotion11 * dt;

#Lagrange Coefficients
var rvar11 = a11 + (r11-a11)*cos(DE[11]) + sigma11*sqrt(a11)*sin(DE[11]);
var F11 = 1 - a11/r11 * (1-cos(DE[11]));
var G11 = a11*sigma11*(1-cos(DE[11])) + r11*sqrt(a11)*sin(DE[11]);
var Ft11 = -sqrt(a11)/(r11*rvar11) * sin(DE[11]);
var Gt11 = 1 - a11/rvar11*(1-cos(DE[11]));

subject to KeplerEquations11: 
	DM11 - DE[11] - sigma11/sqrt(a11) * (1 - cos(DE[11])) + (1 - r11/a11)*sin(DE[11]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 12
var x12 = F11*x11 + G11*(dx11);
var y12 = F11*y11 + G11*(dy11);
var z12 = F11*z11 + G11*(dz11);
var dx12 = Ft11*x11 + Gt11*(dx11)+ux[12];
var dy12 = Ft11*y11 + Gt11*(dy11)+uy[12];
var dz12 = Ft11*z11 + Gt11*(dz11)+uz[12];

#Basic definitions
var r12 = sqrt(x12^2+y12^2+z12^2);
var v12 = sqrt(dx12^2+dy12^2+dz12^2);
var a12 = 1 / (2/r12 - v12^2);
var sigma12 = x12*dx12+y12*dy12+z12*dz12;
var meanmotion12 = sqrt(1/a12^3);
var DM12 = meanmotion12 * dt;

#Lagrange Coefficients
var rvar12 = a12 + (r12-a12)*cos(DE[12]) + sigma12*sqrt(a12)*sin(DE[12]);
var F12 = 1 - a12/r12 * (1-cos(DE[12]));
var G12 = a12*sigma12*(1-cos(DE[12])) + r12*sqrt(a12)*sin(DE[12]);
var Ft12 = -sqrt(a12)/(r12*rvar12) * sin(DE[12]);
var Gt12 = 1 - a12/rvar12*(1-cos(DE[12]));

subject to KeplerEquations12: 
	DM12 - DE[12] - sigma12/sqrt(a12) * (1 - cos(DE[12])) + (1 - r12/a12)*sin(DE[12]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 13
var x13 = F12*x12 + G12*(dx12);
var y13 = F12*y12 + G12*(dy12);
var z13 = F12*z12 + G12*(dz12);
var dx13 = Ft12*x12 + Gt12*(dx12)+ux[13];
var dy13 = Ft12*y12 + Gt12*(dy12)+uy[13];
var dz13 = Ft12*z12 + Gt12*(dz12)+uz[13];

#Basic definitions
var r13 = sqrt(x13^2+y13^2+z13^2);
var v13 = sqrt(dx13^2+dy13^2+dz13^2);
var a13 = 1 / (2/r13 - v13^2);
var sigma13 = x13*dx13+y13*dy13+z13*dz13;
var meanmotion13 = sqrt(1/a13^3);
var DM13 = meanmotion13 * dt;

#Lagrange Coefficients
var rvar13 = a13 + (r13-a13)*cos(DE[13]) + sigma13*sqrt(a13)*sin(DE[13]);
var F13 = 1 - a13/r13 * (1-cos(DE[13]));
var G13 = a13*sigma13*(1-cos(DE[13])) + r13*sqrt(a13)*sin(DE[13]);
var Ft13 = -sqrt(a13)/(r13*rvar13) * sin(DE[13]);
var Gt13 = 1 - a13/rvar13*(1-cos(DE[13]));

subject to KeplerEquations13: 
	DM13 - DE[13] - sigma13/sqrt(a13) * (1 - cos(DE[13])) + (1 - r13/a13)*sin(DE[13]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 14
var x14 = F13*x13 + G13*(dx13);
var y14 = F13*y13 + G13*(dy13);
var z14 = F13*z13 + G13*(dz13);
var dx14 = Ft13*x13 + Gt13*(dx13)+ux[14];
var dy14 = Ft13*y13 + Gt13*(dy13)+uy[14];
var dz14 = Ft13*z13 + Gt13*(dz13)+uz[14];

#Basic definitions
var r14 = sqrt(x14^2+y14^2+z14^2);
var v14 = sqrt(dx14^2+dy14^2+dz14^2);
var a14 = 1 / (2/r14 - v14^2);
var sigma14 = x14*dx14+y14*dy14+z14*dz14;
var meanmotion14 = sqrt(1/a14^3);
var DM14 = meanmotion14 * dt;

#Lagrange Coefficients
var rvar14 = a14 + (r14-a14)*cos(DE[14]) + sigma14*sqrt(a14)*sin(DE[14]);
var F14 = 1 - a14/r14 * (1-cos(DE[14]));
var G14 = a14*sigma14*(1-cos(DE[14])) + r14*sqrt(a14)*sin(DE[14]);
var Ft14 = -sqrt(a14)/(r14*rvar14) * sin(DE[14]);
var Gt14 = 1 - a14/rvar14*(1-cos(DE[14]));

subject to KeplerEquations14: 
	DM14 - DE[14] - sigma14/sqrt(a14) * (1 - cos(DE[14])) + (1 - r14/a14)*sin(DE[14]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 15
var x15 = F14*x14 + G14*(dx14);
var y15 = F14*y14 + G14*(dy14);
var z15 = F14*z14 + G14*(dz14);
var dx15 = Ft14*x14 + Gt14*(dx14)+ux[15];
var dy15 = Ft14*y14 + Gt14*(dy14)+uy[15];
var dz15 = Ft14*z14 + Gt14*(dz14)+uz[15];

#Basic definitions
var r15 = sqrt(x15^2+y15^2+z15^2);
var v15 = sqrt(dx15^2+dy15^2+dz15^2);
var a15 = 1 / (2/r15 - v15^2);
var sigma15 = x15*dx15+y15*dy15+z15*dz15;
var meanmotion15 = sqrt(1/a15^3);
var DM15 = meanmotion15 * dt;

#Lagrange Coefficients
var rvar15 = a15 + (r15-a15)*cos(DE[15]) + sigma15*sqrt(a15)*sin(DE[15]);
var F15 = 1 - a15/r15 * (1-cos(DE[15]));
var G15 = a15*sigma15*(1-cos(DE[15])) + r15*sqrt(a15)*sin(DE[15]);
var Ft15 = -sqrt(a15)/(r15*rvar15) * sin(DE[15]);
var Gt15 = 1 - a15/rvar15*(1-cos(DE[15]));

subject to KeplerEquations15: 
	DM15 - DE[15] - sigma15/sqrt(a15) * (1 - cos(DE[15])) + (1 - r15/a15)*sin(DE[15]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 16
var x16 = F15*x15 + G15*(dx15);
var y16 = F15*y15 + G15*(dy15);
var z16 = F15*z15 + G15*(dz15);
var dx16 = Ft15*x15 + Gt15*(dx15)+ux[16];
var dy16 = Ft15*y15 + Gt15*(dy15)+uy[16];
var dz16 = Ft15*z15 + Gt15*(dz15)+uz[16];

#Basic definitions
var r16 = sqrt(x16^2+y16^2+z16^2);
var v16 = sqrt(dx16^2+dy16^2+dz16^2);
var a16 = 1 / (2/r16 - v16^2);
var sigma16 = x16*dx16+y16*dy16+z16*dz16;
var meanmotion16 = sqrt(1/a16^3);
var DM16 = meanmotion16 * dt;

#Lagrange Coefficients
var rvar16 = a16 + (r16-a16)*cos(DE[16]) + sigma16*sqrt(a16)*sin(DE[16]);
var F16 = 1 - a16/r16 * (1-cos(DE[16]));
var G16 = a16*sigma16*(1-cos(DE[16])) + r16*sqrt(a16)*sin(DE[16]);
var Ft16 = -sqrt(a16)/(r16*rvar16) * sin(DE[16]);
var Gt16 = 1 - a16/rvar16*(1-cos(DE[16]));

subject to KeplerEquations16: 
	DM16 - DE[16] - sigma16/sqrt(a16) * (1 - cos(DE[16])) + (1 - r16/a16)*sin(DE[16]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 17
var x17 = F16*x16 + G16*(dx16);
var y17 = F16*y16 + G16*(dy16);
var z17 = F16*z16 + G16*(dz16);
var dx17 = Ft16*x16 + Gt16*(dx16)+ux[17];
var dy17 = Ft16*y16 + Gt16*(dy16)+uy[17];
var dz17 = Ft16*z16 + Gt16*(dz16)+uz[17];

#Basic definitions
var r17 = sqrt(x17^2+y17^2+z17^2);
var v17 = sqrt(dx17^2+dy17^2+dz17^2);
var a17 = 1 / (2/r17 - v17^2);
var sigma17 = x17*dx17+y17*dy17+z17*dz17;
var meanmotion17 = sqrt(1/a17^3);
var DM17 = meanmotion17 * dt;

#Lagrange Coefficients
var rvar17 = a17 + (r17-a17)*cos(DE[17]) + sigma17*sqrt(a17)*sin(DE[17]);
var F17 = 1 - a17/r17 * (1-cos(DE[17]));
var G17 = a17*sigma17*(1-cos(DE[17])) + r17*sqrt(a17)*sin(DE[17]);
var Ft17 = -sqrt(a17)/(r17*rvar17) * sin(DE[17]);
var Gt17 = 1 - a17/rvar17*(1-cos(DE[17]));

subject to KeplerEquations17: 
	DM17 - DE[17] - sigma17/sqrt(a17) * (1 - cos(DE[17])) + (1 - r17/a17)*sin(DE[17]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 18
var x18 = F17*x17 + G17*(dx17);
var y18 = F17*y17 + G17*(dy17);
var z18 = F17*z17 + G17*(dz17);
var dx18 = Ft17*x17 + Gt17*(dx17)+ux[18];
var dy18 = Ft17*y17 + Gt17*(dy17)+uy[18];
var dz18 = Ft17*z17 + Gt17*(dz17)+uz[18];

#Basic definitions
var r18 = sqrt(x18^2+y18^2+z18^2);
var v18 = sqrt(dx18^2+dy18^2+dz18^2);
var a18 = 1 / (2/r18 - v18^2);
var sigma18 = x18*dx18+y18*dy18+z18*dz18;
var meanmotion18 = sqrt(1/a18^3);
var DM18 = meanmotion18 * dt;

#Lagrange Coefficients
var rvar18 = a18 + (r18-a18)*cos(DE[18]) + sigma18*sqrt(a18)*sin(DE[18]);
var F18 = 1 - a18/r18 * (1-cos(DE[18]));
var G18 = a18*sigma18*(1-cos(DE[18])) + r18*sqrt(a18)*sin(DE[18]);
var Ft18 = -sqrt(a18)/(r18*rvar18) * sin(DE[18]);
var Gt18 = 1 - a18/rvar18*(1-cos(DE[18]));

subject to KeplerEquations18: 
	DM18 - DE[18] - sigma18/sqrt(a18) * (1 - cos(DE[18])) + (1 - r18/a18)*sin(DE[18]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 19
var x19 = F18*x18 + G18*(dx18);
var y19 = F18*y18 + G18*(dy18);
var z19 = F18*z18 + G18*(dz18);
var dx19 = Ft18*x18 + Gt18*(dx18)+ux[19];
var dy19 = Ft18*y18 + Gt18*(dy18)+uy[19];
var dz19 = Ft18*z18 + Gt18*(dz18)+uz[19];

#Basic definitions
var r19 = sqrt(x19^2+y19^2+z19^2);
var v19 = sqrt(dx19^2+dy19^2+dz19^2);
var a19 = 1 / (2/r19 - v19^2);
var sigma19 = x19*dx19+y19*dy19+z19*dz19;
var meanmotion19 = sqrt(1/a19^3);
var DM19 = meanmotion19 * dt;

#Lagrange Coefficients
var rvar19 = a19 + (r19-a19)*cos(DE[19]) + sigma19*sqrt(a19)*sin(DE[19]);
var F19 = 1 - a19/r19 * (1-cos(DE[19]));
var G19 = a19*sigma19*(1-cos(DE[19])) + r19*sqrt(a19)*sin(DE[19]);
var Ft19 = -sqrt(a19)/(r19*rvar19) * sin(DE[19]);
var Gt19 = 1 - a19/rvar19*(1-cos(DE[19]));

subject to KeplerEquations19: 
	DM19 - DE[19] - sigma19/sqrt(a19) * (1 - cos(DE[19])) + (1 - r19/a19)*sin(DE[19]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 20
var x20 = F19*x19 + G19*(dx19);
var y20 = F19*y19 + G19*(dy19);
var z20 = F19*z19 + G19*(dz19);
var dx20 = Ft19*x19 + Gt19*(dx19)+ux[20];
var dy20 = Ft19*y19 + Gt19*(dy19)+uy[20];
var dz20 = Ft19*z19 + Gt19*(dz19)+uz[20];

#Basic definitions
var r20 = sqrt(x20^2+y20^2+z20^2);
var v20 = sqrt(dx20^2+dy20^2+dz20^2);
var a20 = 1 / (2/r20 - v20^2);
var sigma20 = x20*dx20+y20*dy20+z20*dz20;
var meanmotion20 = sqrt(1/a20^3);
var DM20 = meanmotion20 * dt;

#Lagrange Coefficients
var rvar20 = a20 + (r20-a20)*cos(DE[20]) + sigma20*sqrt(a20)*sin(DE[20]);
var F20 = 1 - a20/r20 * (1-cos(DE[20]));
var G20 = a20*sigma20*(1-cos(DE[20])) + r20*sqrt(a20)*sin(DE[20]);
var Ft20 = -sqrt(a20)/(r20*rvar20) * sin(DE[20]);
var Gt20 = 1 - a20/rvar20*(1-cos(DE[20]));

subject to KeplerEquations20: 
	DM20 - DE[20] - sigma20/sqrt(a20) * (1 - cos(DE[20])) + (1 - r20/a20)*sin(DE[20]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 21
var x21 = F20*x20 + G20*(dx20);
var y21 = F20*y20 + G20*(dy20);
var z21 = F20*z20 + G20*(dz20);
var dx21 = Ft20*x20 + Gt20*(dx20)+ux[21];
var dy21 = Ft20*y20 + Gt20*(dy20)+uy[21];
var dz21 = Ft20*z20 + Gt20*(dz20)+uz[21];

#Basic definitions
var r21 = sqrt(x21^2+y21^2+z21^2);
var v21 = sqrt(dx21^2+dy21^2+dz21^2);
var a21 = 1 / (2/r21 - v21^2);
var sigma21 = x21*dx21+y21*dy21+z21*dz21;
var meanmotion21 = sqrt(1/a21^3);
var DM21 = meanmotion21 * dt;

#Lagrange Coefficients
var rvar21 = a21 + (r21-a21)*cos(DE[21]) + sigma21*sqrt(a21)*sin(DE[21]);
var F21 = 1 - a21/r21 * (1-cos(DE[21]));
var G21 = a21*sigma21*(1-cos(DE[21])) + r21*sqrt(a21)*sin(DE[21]);
var Ft21 = -sqrt(a21)/(r21*rvar21) * sin(DE[21]);
var Gt21 = 1 - a21/rvar21*(1-cos(DE[21]));

subject to KeplerEquations21: 
	DM21 - DE[21] - sigma21/sqrt(a21) * (1 - cos(DE[21])) + (1 - r21/a21)*sin(DE[21]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 22
var x22 = F21*x21 + G21*(dx21);
var y22 = F21*y21 + G21*(dy21);
var z22 = F21*z21 + G21*(dz21);
var dx22 = Ft21*x21 + Gt21*(dx21)+ux[22];
var dy22 = Ft21*y21 + Gt21*(dy21)+uy[22];
var dz22 = Ft21*z21 + Gt21*(dz21)+uz[22];

#Basic definitions
var r22 = sqrt(x22^2+y22^2+z22^2);
var v22 = sqrt(dx22^2+dy22^2+dz22^2);
var a22 = 1 / (2/r22 - v22^2);
var sigma22 = x22*dx22+y22*dy22+z22*dz22;
var meanmotion22 = sqrt(1/a22^3);
var DM22 = meanmotion22 * dt;

#Lagrange Coefficients
var rvar22 = a22 + (r22-a22)*cos(DE[22]) + sigma22*sqrt(a22)*sin(DE[22]);
var F22 = 1 - a22/r22 * (1-cos(DE[22]));
var G22 = a22*sigma22*(1-cos(DE[22])) + r22*sqrt(a22)*sin(DE[22]);
var Ft22 = -sqrt(a22)/(r22*rvar22) * sin(DE[22]);
var Gt22 = 1 - a22/rvar22*(1-cos(DE[22]));

subject to KeplerEquations22: 
	DM22 - DE[22] - sigma22/sqrt(a22) * (1 - cos(DE[22])) + (1 - r22/a22)*sin(DE[22]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 23
var x23 = F22*x22 + G22*(dx22);
var y23 = F22*y22 + G22*(dy22);
var z23 = F22*z22 + G22*(dz22);
var dx23 = Ft22*x22 + Gt22*(dx22)+ux[23];
var dy23 = Ft22*y22 + Gt22*(dy22)+uy[23];
var dz23 = Ft22*z22 + Gt22*(dz22)+uz[23];

#Basic definitions
var r23 = sqrt(x23^2+y23^2+z23^2);
var v23 = sqrt(dx23^2+dy23^2+dz23^2);
var a23 = 1 / (2/r23 - v23^2);
var sigma23 = x23*dx23+y23*dy23+z23*dz23;
var meanmotion23 = sqrt(1/a23^3);
var DM23 = meanmotion23 * dt;

#Lagrange Coefficients
var rvar23 = a23 + (r23-a23)*cos(DE[23]) + sigma23*sqrt(a23)*sin(DE[23]);
var F23 = 1 - a23/r23 * (1-cos(DE[23]));
var G23 = a23*sigma23*(1-cos(DE[23])) + r23*sqrt(a23)*sin(DE[23]);
var Ft23 = -sqrt(a23)/(r23*rvar23) * sin(DE[23]);
var Gt23 = 1 - a23/rvar23*(1-cos(DE[23]));

subject to KeplerEquations23: 
	DM23 - DE[23] - sigma23/sqrt(a23) * (1 - cos(DE[23])) + (1 - r23/a23)*sin(DE[23]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 24
var x24 = F23*x23 + G23*(dx23);
var y24 = F23*y23 + G23*(dy23);
var z24 = F23*z23 + G23*(dz23);
var dx24 = Ft23*x23 + Gt23*(dx23)+ux[24];
var dy24 = Ft23*y23 + Gt23*(dy23)+uy[24];
var dz24 = Ft23*z23 + Gt23*(dz23)+uz[24];

#Basic definitions
var r24 = sqrt(x24^2+y24^2+z24^2);
var v24 = sqrt(dx24^2+dy24^2+dz24^2);
var a24 = 1 / (2/r24 - v24^2);
var sigma24 = x24*dx24+y24*dy24+z24*dz24;
var meanmotion24 = sqrt(1/a24^3);
var DM24 = meanmotion24 * dt;

#Lagrange Coefficients
var rvar24 = a24 + (r24-a24)*cos(DE[24]) + sigma24*sqrt(a24)*sin(DE[24]);
var F24 = 1 - a24/r24 * (1-cos(DE[24]));
var G24 = a24*sigma24*(1-cos(DE[24])) + r24*sqrt(a24)*sin(DE[24]);
var Ft24 = -sqrt(a24)/(r24*rvar24) * sin(DE[24]);
var Gt24 = 1 - a24/rvar24*(1-cos(DE[24]));

subject to KeplerEquations24: 
	DM24 - DE[24] - sigma24/sqrt(a24) * (1 - cos(DE[24])) + (1 - r24/a24)*sin(DE[24]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node 25
var x25 = F24*x24 + G24*(dx24);
var y25 = F24*y24 + G24*(dy24);
var z25 = F24*z24 + G24*(dz24);
var dx25 = Ft24*x24 + Gt24*(dx24)+ux[25];
var dy25 = Ft24*y24 + Gt24*(dy24)+uy[25];
var dz25 = Ft24*z24 + Gt24*(dz24)+uz[25];

#Basic definitions
var r25 = sqrt(x25^2+y25^2+z25^2);
var v25 = sqrt(dx25^2+dy25^2+dz25^2);
var a25 = 1 / (2/r25 - v25^2);
var sigma25 = x25*dx25+y25*dy25+z25*dz25;
var meanmotion25 = sqrt(1/a25^3);
var DM25 = meanmotion25 * dt;

#Lagrange Coefficients
var rvar25 = a25 + (r25-a25)*cos(DE[25]) + sigma25*sqrt(a25)*sin(DE[25]);
var F25 = 1 - a25/r25 * (1-cos(DE[25]));
var G25 = a25*sigma25*(1-cos(DE[25])) + r25*sqrt(a25)*sin(DE[25]);
var Ft25 = -sqrt(a25)/(r25*rvar25) * sin(DE[25]);
var Gt25 = 1 - a25/rvar25*(1-cos(DE[25]));

subject to KeplerEquations25: 
	DM25 - DE[25] - sigma25/sqrt(a25) * (1 - cos(DE[25])) + (1 - r25/a25)*sin(DE[25]) = 0;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Node n: Arrival node
var xn = F25*x25 + G25*dx25;
var yn = F25*y25 + G25*dy25;
var zn = F25*z25 + G25*dz25;
var dxn = Ft25*x25 + Gt25*dx25 + ux[25];
var dyn = Ft25*y25 + Gt25*dy25 + uy[25];
var dzn = Ft25*z25 + Gt25*dz25 + uz[25];

#Basic definitions
var rn = sqrt(xn^2+yn^2+zn^2);
var vn = sqrt(dxn^2+dyn^2+dzn^2);
var an = 1 / (2/rn - vn^2);
var sigman = xn*dxn+yn*dyn+zn*dzn;
var meanmotionn = sqrt(1/an^3);
var DMn = meanmotionn * dt;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Match Constraint
subject to 
	FinalPositionx  : xn = xf;
	FinalPositiony  : yn = yf;
	FinalPositionz  : zn = zf;
	FinalVelocityx  : dxn = dxf - ux[n];
	FinalVelocityy  : dyn = dyf - uy[n];
	FinalVelocityz  : dzn = dzf - uz[n];
	
#--------------------------------------------------------------------------
