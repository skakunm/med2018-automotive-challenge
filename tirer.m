function [Fx, Fy] = tirer(k, a, Fz)
Fz0 = 3252;
PCX1=1.63;
PDX1=1.06;
PDX2=-0.0492;
PEX1=0.5;
PEX2=-0.11;
PEX3=-0.06;
PEX4=0;
PKX1=19.7;
PKX2=-0.15;
PKX3=0.18;
PHX1=-0.0005;
PHX2=8.5e-5;
PVX1=0;
PVX2=0;
RBX1=9.0;
RBX2=-8.6;
RCX1=1.131;
REX1=0.081;
REX2=-0.15;
RHX1=-0.029;
PCY1=1.28;
PDY1=-0.92;
PDY2=0.22;
PEY1=-1.1;
PEY2=0.65;
PEY3=-0.65;
PKY1=-13.06;
PKY2=1.77;
PHY1=0.0034;
PHY2=-0.003;
PVY1=0.044;
PVY2=-0.030;
RBY1=6.4;
RBY2=7.91;
RBY3=-0.059;
RCY1=1.16;
REY1=0.22;
REY2=0.43;
RHY1=0.0007;
RHY2=0.023;
RVY1=0;
RVY2=0;
RVY4=10;
RVY5=1.94;
RVY6=-50;
dfz = (Fz - Fz0)/Fz0;

Shx = PHX1 + PHX2*dfz;
Svx = Fz * (PVX1 + PVX2 * dfz);
kx = k + Shx;
Cx = PCX1;
Dx = Fz*(PDX1 + PDX2 * dfz);
Ex = (PEX1 + PEX2 * dfz + PEX3 * dfz^2) * (1 - PEX4 * sign(kx));
Bx = Fz * (PKX1 + PKX2 * dfz) * exp(PKX3 * dfz) / (Cx*Dx);
Fx0 = Dx*sin(Cx*atan(Bx*kx - Ex*(Bx*kx - atan(Bx*kx)))) + Svx;

Bxa = RBX1 * cos(atan(RBX2 * k));
Cxa = RCX1;
Exa = REX1 + REX2*dfz;
Shxa = RHX1;
as = a + Shxa;
Gxa = (cos(Cxa*atan(Bxa*as - Exa*(Bxa*as - atan(Bxa*as)))))/...
(cos(Cxa*atan(Bxa*Shxa - Exa*(Bxa*Shxa - atan(Bxa*Shxa)))));

Fx = Fx0 * Gxa;


Shy = PHY1 + PHY2*dfz;
Svy = Fz * (PVY1 + PVY2 * dfz);
ay = a + Shy;
Cy = PCY1;
Dy = Fz*(PDY1 + PDY2 * dfz);
Ey = (PEY1 + PEY2 * dfz) * (1 - PEY3 * sign(ay));
By = PKY1 * Fz0 * sin(2 * atan(Fz/(PKY2 * Fz0)))  / (Cy*Dy);
Fy0 = Dy*sin(Cy*atan(By*ay - Ey*(By*ay - atan(By*ay)))) + Svy;

Svyk = (PDY1 + PDY2*dfz)*Fz*(RVY1 + RVY2*dfz)*...
    cos(atan(RVY4*a))*sin(RVY5 * atan(RVY6*k));
Shyk = RHY1 + RHY2 * dfz;
ks = k + Shyk;
Byk = RBY1 * cos(atan(RBY2 * (a - RBY3)));
Cyk = RCY1;
Eyk = REY1 + REY2*dfz;
Gyk = (cos(Cyk*atan(Byk*ks - Eyk*(Byk*ks - atan(Byk*ks)))))/...
(cos(Cyk*atan(Byk*Shyk - Eyk*(Byk*Shyk - atan(Byk*Shyk)))));

Fy = Fy0 * Gyk + Svyk;
end