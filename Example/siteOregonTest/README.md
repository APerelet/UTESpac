# UTESpac.m Settings

```matlab
%path to folder containing the site* folder
info.rootFolder = '/UTESpac/Example/';

info.foldStruct = '';

info.FileForm = '(?<serial>\d+)[.](?<TableName>\w*)_(?<Year>\d{4})_(?<Month>\d{2})_(?<Day>\d{2})_(?<Hour>\d{2})(?<Minute>\d{2}).dat';

%Use coefficients from PFInfo.m
%If PFInfo.m does not exist, this should be local
info.PF.globalCalculation = 'local'; 
                                      
%if PFInfo.m does not exist and global planar fit is desired this should be true
info.PF.recalculateGlobalCoefficients = false;  

info.PF.avgPer = 30;
```
Check Template to make sure it matches Header files - case sensitive

```
template.u = 'Ux_*'; % sonic u  --   [m/s]
template.v = 'Uy_*'; % sonic v  --   [m/s]
template.w = 'Uz_*'; % sonic w  --   [m/s]
template.Tson = 'T_Sonic_*'; % sonic T  --   [C or K]
template.sonDiagnostic = 'sonic_diag_*'; % sonic diagnostic  --  [-]
template.fw = 'fw_*'; % sonic finewires to be used for Eddy Covariance  --  [C]
template.RH = 'HMP_RH_*'; % slow response relative humidity for virtual temperature calculation  --  [Fract or %]
template.T = 'HMP_T_*'; % slow response temperature  --  [C]
template.P = 'Pressure'; % pressure  --  [kPa or mBar]
template.irgaH2O = 'H2O_*'; % for use with Campbell EC150 and IRGASON.  WPL corrections applied  --  [g/m^3]
template.irgaH2OsigStrength = 'H2Osig_*'; % EC150 Signal Strength  --  [-]
template.irgaCO2 = 'CO2_*'; % for use with Campbell EC150 and IRGASON.  WPL corrections applied  --  [mg/m^3]
template.irgaCO2sigStrength = 'CO2sig_*'; % EC150 Signal Strength  --  [-]
template.irgaGasDiag = 'gas_diag_*'; % EC150 gas (CO2 and H2O) diagnostic, 0-> Okay  --  [-]
template.LiH2O = 'LiH2O_*'; % for use with Licor 7500.  WPL corrections applied  --  [mmol/mol]
template.LiCO2 = 'LiCO2_*'; % for use with Licor 7500.  WPL corrections applied  --  [mmol/mol]
template.LiGasDiag = 'Li_gas_diag_*'; % Li7500 gas (CO2 and H2O) diagnostic >~230 -> Okay  --  [-]
template.KH2O = 'KH2O_H2O_*'; % for use with KH2Os.  WPL and O2 corrections applied  --  [g/m^3]
template.cup = 'cup_*';  % wind speed from cup anemometers  --  [m/s]
template.birdSpd = 'wbSpd_*';  % wind speed from prop anemometer  --  [m/s]
template.birdDir = 'wbDir_*';  % wind direction from vain or prop anemometer  --  [deg]
```

When running UTESpac.m make sure found instruments contain all the instruments mentioned in header files
