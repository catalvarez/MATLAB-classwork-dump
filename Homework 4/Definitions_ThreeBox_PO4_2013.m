%   This is Definitions_ThreeBox_PO4.m
%   The idea here is to solve the Sarmiento and Toggweiller 1984 3 box
%   model system for phosphate only, using a first order kinetics rule for
%   bio productivity.

%   Give each variable_box # its own number for clarity in ode vector
%   calls

po4_1 = 1; % Low Lat. Surface Box
po4_2 = 2; % High Lat. Surface Box
po4_3 = 3; % Deep Box

%   For running 'dafunPE_jfa.m'

global phflag k1k2flag Ca Mg;

Ca = 10.3e-3;						% 10.3 (mol/kg) modern
Mg = 53.0e-3;						% 53.0 (mol/kg) modern
phflag = 0;							% flag for using the total pH scale in script
k1k2flag = 0;						% flag for using the Meybach acid ks in script

%	Water Fluxes

fLH = 1e6;							% High Lat-Low Lat Exchange.  Value in m^3/sec
fHD = 40e6;							% High Lat-Deep Exchange.  Value in m^3/sec
fLD = 1e6;							% Low Lat-Deep Exchange.  Value in m^3/sec
T = 25e6;								% Overturning.  Value in m^3/sec
SecsPerYear = 3.14e7;				% Seconds in a year.

%	Volumes

AreaOcean = 361e12;                     % Ocean area in m^2
DepthOcean = 3780;                      % Ocean depth in m
VTot = AreaOcean*DepthOcean;	% Ocean volume as a product of area and mean depth in m^3.
AreafracH = 0.15;
AreafracL = 0.85;

VH = AreafracH*AreaOcean*250;   % Volume of the High Lat box, 250m deep
VL = AreafracL*AreaOcean*100;   % Volume of the Low Lat box, 100m deep
VD = VTot-(VH+VL);                      % Volume of the Deep box as a difference from other two
OceanVolArray = [VL VH VD];         	% array of box volumes to use on rhs of dcdt expression
InvOceanVolArray = 1./OceanVolArray;

AgeD_fHD = VD/(fHD)/SecsPerYear;
AgeD_T = VD/(T)/SecsPerYear;

%	Biological Reaction Rate Constants,  Total PO4, and remin depth fraction

kH = 1.0e-9;							% Rate constant for high lat south Biological prod.  1/sec.
kL = 6.0e-5;							% Rate constant for low lat Biological prod.  1/sec.

%PO4ConcMean = 2.1/1000;					% Mean concentration of PO4 in thte ocean (µmole/kg converted to mole/m^3);