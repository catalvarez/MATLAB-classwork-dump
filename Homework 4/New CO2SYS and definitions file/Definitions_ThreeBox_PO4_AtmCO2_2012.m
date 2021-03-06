%   This is Definitions_ThreeBox_PO4.m
%   The idea here is to solve the Sarmiento and Toggweiller 1984 3 box
%   model system for phosphate only, using a first order kinetics rule for
%   bio productivity.

%   Give each variable_box # its own number for clarity in ode vector
%   calls

po4_1 = 1;  %low
po4_2 = 2;  %high
po4_3 = 3;  %deep
po4_4 = 4;  %atmo
alk_1 = 5;
alk_2 = 6;
alk_3 = 7;
alk_4 = 8;
dic_1 = 9;
dic_2 = 10;
dic_3 = 11;
pco2_4 = 12;

%   For running 'dafunPE_jfa.m'

global phflag k1k2flag Ca Mg;

Ca = 10.3e-3;						% 10.3 (mol/kg) modern
Mg = 53.0e-3;						% 53.0 (mol/kg) modern
phflag = 0;							% flag for using the total pH scale in script
k1k2flag = 0;						% flag for using the Meybach acid ks in script

%	Water Fluxes

fLH = 1e6;							% High Lat-Low Lat Exchange.  Value in m^3/sec
fHD = 3e6;							% High Lat-Deep Exchange.  Value in m^3/sec
fLD = 1e6;							% Low Lat-Deep Exchange.  Value in m^3/sec
T = 25e6;								% Overturning.  Value in m^3/sec
SecsPerYear = 3.14e7;				% Seconds in a year.

%	Volumes

MassAtm = 5.148e18;                 % Atm mass in kg
DensAtm = 1.2041;                   % Dens atm at sea level (kg/m^3)
V4 = MassAtm./DensAtm;				% Volume of the atmosphere (m^3)
AreaAtm = 5.10072e8;                % Area of atm in km^2
AreaAtm = AreaAtm*1000000;          % Area of atm in m^2
R = 0.082053;                       % Gas constant (l atm/�K/mol)
R = R./1000;                        % Gas constant (m^3 atm/�K/mol)
Tatm = 25+273.15;                   % Temperature atm

AreaOcean = 361e12;                     % Ocean area in m^2
DepthOcean = 3780;                      % Ocean depth in m
VTot = AreaOcean*DepthOcean;	% Ocean volume as a product of area and mean depth in m^3.
AreafracH = 0.15;
AreafracL = 0.85;

Area_1 = AreafracL.*AreaOcean;
Area_2 = AreafracH.*AreaOcean;
VH = AreafracH*AreaOcean*250;   % Volume of the High Lat box, 250m deep
VL = AreafracL*AreaOcean*100;   % Volume of the Low Lat box, 100m deep
VD = VTot-(VH+VL);                      % Volume of the Deep box as a difference from other two
OceanVolArray = [VL VH VD V4/(R*Tatm)];         	% array of box volumes to use on rhs of dcdt expression
InvOceanVolArray = 1./OceanVolArray;

%AgeD_fHD = VD/(fHD)/SecsPerYear;
%AgeD_T = VD/(T)/SecsPerYear;

%	Temperatures and pressures of boxes

TC = [25 6 2]';				% Temperature variable for Zeebe script
P = [100/20 250/20 350]';		% Pressure variable for Zeebe script


%	Biological Reaction Rate Constants,  Total PO4, and remin depth fraction

kH = 1.0e-9;							% Rate constant for high lat south Biological prod.  1/sec.
kL = 6.0e-5;							% Rate constant for low lat Biological prod.  1/sec.

%PO4ConcMean = 2.1/1000;					% Mean concentration of PO4 in thte ocean (�mole/kg converted to mole/m^3);

%	Redfield Ratios

rdic_po4 = 131;						% Redfiled ratio with 106 organic C and 25 CaCO3 C
ralk_po4 = 34;						% Redfield ratio from Sarmiento and Toggweiler

%   Gas Exchange Stuff

PV_1 = 3/24/60/60;                  % Piston velocity box #1 3m/day
PV_2 = 3/24/60/60;                  % Piston velocity box #2 3m/day

B_1 = 30;                         % Solubility (beta) box #1 (mol/m^3/atm)
B_2 = 60;                         % Solubility (beta) box #2 (mol/m^3/atm)
