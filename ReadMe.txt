
================================================================================
CET86: Complex Chemical Equilibrium and Thermal Transport Program (WINDOWS-PORT)
================================================================================

Program Capabilities:

Thermodynamic and Thermal Transport Mixture Properties:

Chemical equilibrium compositions are obtained by the method of free energy minimization. A thermodynamic state is characterized by two independent state variables, such as temperature and pressure. If pressure is 'one of the state variables, Gibbs energy is minimized. This is the case for the following combination of variables permitted to be assigned by the program: temperature and pressure (tp), enthalpy and pressure (hp), and entropy and pressure (sp). If volume (or density) is one of the state variables, Helmholtz energy is minimized. This is the case for the following combination of variables permitted by the program: temperature and volume (or density) (tv), internal energy and volume (or density) (uv), and entropy and volume (or density) (sv).

It is assumed that all gases are ideal and that interactions among phases can be neglected. An ideal-gas equation of state is used to represent the mixture and is assumed to be correct even when small amounts of condensed phases are present. Equilibrium properties of plasmas (mixtures containing ionized species) may also be calculated if the plasma is considered to be ideal, that is, if columbic interactions are not considered.

Thermodynamic properties of mixtures include the contribution of condensed as well as gaseous phases. However, thermal transport mixture properties include the contributions of gas-phase species only. If condensed phases are present, mole fractions for the gas-phase species are first normalized to gases only prior to calculating thermal transport mixture properties. The thermodynamic and thermal transport mixture properties calculated by the program are discussed in the section


1. Pre-Requisits/Reference: (https://ntrs.nasa.gov/search.jsp?R=19940028442)

2. Documentation: (https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19940028442.pdf)

3. Example problem located in "./examples/":

4. Example problem: "Combustion Temperatures for Nitromethane & Nitrous Oxide at various pressures and Oxidizer-to-Fuel Ratios" located: "./examples/"
4a. INPUT FILE #1: "Nitromethane_NitrousOxide.inp".

==================================================
START: Input File: "Nitromethane_NitrousOxide.inp"
==================================================
THERMO                                                                          
   300.000  1000.000  5000.000                                                  
E                 L02/67E  1.0  0.0  0.0  0.G   300.000  5000.000      0.00000 1
 0.25000000E 01     0.00000000     0.00000000     0.00000000     0.00000000    2
-0.74537496E 03-0.11734026E 02 0.25000000E 01     0.00000000     0.00000000    3
     0.00000000     0.00000000-0.74537500E 03-0.11734026E 02     0.00000000    4
CH3NO             J10/82C  1.H  3.N  1.O  1.G   300.000  5000.000              1
 0.41559029E 01 0.11491477E-01-0.45591541E-05 0.82123708E-09-0.55272623E-13    2
 0.61153828E 04 0.47871399E 01 0.28445129E 01 0.13019577E-01-0.32880289E-05    3
-0.14029948E-08 0.69745403E-12 0.66445273E 04 0.12235489E 02                   4
CH3ONO            J10/82C  1.H  3.O  2.N  1.G   300.000  5000.000              1
 0.64071016E 01 0.12673147E-01-0.51597244E-05 0.94678620E-09-0.64593404E-13    2
-0.10759027E 05-0.69541931E 01 0.22027159E 01 0.20894911E-01-0.75951830E-05    3
-0.23690805E-08 0.17073842E-11-0.93791328E 04 0.15746949E 02                   4

... <content removed> ...

ZRO2(B)           J12/65ZR 1.O  2.   0.   0.S  1478.000  2950.000    123.21880 1
 0.89573629E 01     0.00000000     0.00000000     0.00000000     0.00000000    2
-0.13414354E 06-0.45274017E 02 0.89573629E 01     0.00000000     0.00000000    3
     0.00000000     0.00000000-0.13414354E 06-0.45274017E 02     0.00000000    4
ZRO2(L)           J12/65ZR 1.O  2.   0.   0.L  2950.000  5000.000    123.21880 1
 0.10567675E 02     0.00000000     0.00000000     0.00000000     0.00000000    2
-0.12842745E 06-0.54592264E 02 0.10567675E 02     0.00000000     0.00000000    3
     0.00000000     0.00000000-0.12842745E 06-0.54592264E 02     0.00000000    4
END                                                                             
TRANSPORT PROPERTY COEFFICIENTS                                                 
AL                                V2C2  GORDON; NASA TM 86885, OCT 1984         
 V  300.000 1000.000 0.10752557E 01 0.19889058E 03-0.12117144E 05-0.21520631E 01
 V 1000.000 5000.000 0.71350606E 00-0.11856849E 04 0.54275069E 06 0.11828645E 01
 C  300.000 1000.000 0.10752525E 01 0.19888814E 03-0.12116940E 05-0.20074452E 01
 C 1000.000 5000.000 0.71350537E 00-0.11856885E 04 0.54275195E 06 0.13274647E 01
ALCL                              V2C2  GORDON; NASA TM86885, OCT 1984          
 V  300.000 1000.000 0.10793661E 01 0.29479492E 02 0.34836606E 04-0.16604981E 01
 V 1000.000 5000.000 0.56571504E 00-0.61915065E 03 0.84747061E 05 0.24526497E 01
 C  300.000 1000.000 0.98944147E 00-0.77767293E 02 0.95232979E 04-0.10951633E 01
 C 1000.000 5000.000 0.92919002E 00 0.55439951E 03-0.38427598E 06-0.92595436E 00

... <content removed> ...

ZN                                V2C2  GORDON; NASA TM86885, OCT 1984          
 V  300.000 1000.000 0.12002271E 01 0.17586903E 03-0.83995383E 04-0.19217073E 01
 V 1000.000 5000.000 0.53031143E 00-0.11169229E 04 0.30819706E 06 0.36814803E 01
 C  300.000 1000.000 0.12002288E 01 0.17587033E 03-0.83996394E 04-0.26621911E 01
 C 1000.000 5000.000 0.53031130E 00-0.11169233E 04 0.30819716E 06 0.29410104E 01
HE              AR                V1C0                                          
 V  300.000 5000.000 0.47903400E 00-0.24133330E 03 0.34125770E 05 0.27830000E 01
AR              KR                V1C0                                          
 V  300.000 5000.000 0.53955200E 00-0.14537710E 03 0.77105300E 04 0.27820000E 01
CH4             CF4               V1C0                                          
 V  300.000 5000.000 0.13074500E 00-0.55907700E 03 0.55942230E 05 0.52550000E 01
LAST
REACTANTS                                                                       
C 1.     H 3.     N 1.     O 2.     00       100.             G298.15  F
N 2.     O 1.                       00       100.             G298.15  O

NAMELISTS                                                                       
 &INPT2  KASE=123,HP=T,OF=T,MIX=2.5,2.0,1.5,1.0,0.9,0.8,0.7,0.6,0.5,P=100,10,1,TRNSPT=T,TRPACC=.99999&END       
==================================================
END: Input File: "Nitromethane_NitrousOxide.inp"
==================================================


5. Operation: Run "./bin/CET86.exe"


6. OUTPUT FILE: "cet86.out".


==================================
START: Output File: "cet86.out"
==================================
 THERMO                                                                     
 TRANSPORT PR   RTY COEFFICIENTS                                            
 REACTANTS                                                                  
 C   1.0000  H   3.0000  N   1.0000  O   2.0000  00  0.0000  100.000000          0.00  G    298.150     F
 N   2.0000  O   1.0000      0.0000      0.0000  00  0.0000  100.000000          0.00  G    298.150     O
     
    
 NAMELISTS                                                                  
 &INPT2
 KASE    =         123,
 T       = 26*0.0000000E+00  ,
 P       =   100.0000    ,   10.00000    ,   1.000000    , 23*0.0000000E+00  ,
 PSIA    = F,
 MMHG    = F,
 NSQM    = F,
 V       = 26*0.0000000E+00  ,
 RHO     =   100.0000    ,   10.00000    ,   1.000000    , 23*0.0000000E+00  ,
 ERATIO  = F,
 OF      = T,
 FPCT    = F,
 FA      = F,
 MIX     =   2.500000    ,   2.000000    ,   1.500000    ,   1.000000    ,  0.9000000    ,  0.8000000    ,  0.7000000    ,
   0.6000000    ,  0.5000000    , 17*-1.000000      ,
 TP      = F,
 HP      = T,
 SP      = F,
 TV      = F,
 UV      = F,
 SV      = F,
 RKT     = F,
 SHOCK   = F,
 DETN    = F,
 TRACE   =  0.000000000000000E+000,
 S0      =  0.000000000000000E+000,
 SO      =  0.0000000E+00,
 IONS    = F,
 IDEBUG  =           0,
 PHI     = F,
 SIUNIT  = F,
 INHG    = F,
 TRNSPT  = T,
 TRPACC  =  0.999990000000000     ,
 DIF     = F,
 NODATA  = F,
 U       =  1.000000000000000E+030,
 H       =  1.000000000000000E+030
 /
0SPECIES BEING CONSIDERED IN THIS SYSTEM  
  J10/82  CH3NO             J10/82  CH3ONO            J10/82  HCO2              J10/82  CH2NO2            J10/82  CH3NO2          
  J 3/78  C                 J12/67  CH                J12/72  CH2               J 3/61  FORMALDEHYDE      L 4/85  FORMIC ACID     

... <removed content> ...

0NOTE. WEIGHT FRACTION OF FUEL IN TOTAL FUELS AND OF OXIDANT IN TOTAL OXIDANTS
1     
0OF =   0.500000
                            EFFECTIVE FUEL          EFFECTIVE OXIDANT            MIXTURE
 ENTHALPY                       HPP(2)                   HPP(1)                   HSUB0
 (KG-MOL)(DEG K)/KG        -0.12884392E+03           0.22298778E+03          -0.11566685E+02
0KG-FORM.WT./KG               B0P(I,2)                 B0P(I,1)                  B0(I)
 CO                         0.24573969E-01          -0.22720663E-01           0.88090918E-02
 H2O                        0.24573969E-01           0.00000000E+00           0.16382646E-01
 N2                         0.81913232E-02           0.22720663E-01           0.13034437E-01
 CO2                       -0.81913232E-02           0.22720663E-01           0.21126724E-02
   1    4 3366.64  -30.715         -35.170         -25.363         -46.550
   2    4 3152.96  -33.015         -37.791         -27.438         -49.485
   3    4 2921.30  -35.356         -40.517         -29.476         -52.588
1                                    THERMODYNAMIC EQUILIBRIUM COMBUSTION PROPERTIES AT ASSIGNED
0                                                               PRESSURES

 CASE NO.     123
                                                                             WT FRACTION    ENERGY   STATE   TEMP
          CHEMICAL FORMULA                                                   (SEE NOTE)     CAL/MOL          DEG K
 FUEL    C  1.00000   H  3.00000   N  1.00000   O  2.00000                     1.000000   -15628.584   G    298.15
 OXIDANT N  2.00000   O  1.00000                                               1.000000    19502.928   G    298.15
0               O/F=  0.5000    PERCENT FUEL=  66.6667    EQUIVALENCE RATIO= 1.2995    PHI= 2.1631
0THERMODYNAMIC PROPERTIES

 P, ATM            100.00   10.000   1.0000
 T, DEG K          3366.6   3153.0   2921.3
 RHO, G/CC       8.7410-3 9.1633-4 9.6930-5
 H, CAL/G         -22.985  -22.985  -22.985
 U, CAL/G         -300.04  -287.27  -272.83
 G, CAL/G        -8534.62 -8597.21 -8536.73
 S, CAL/(G)(K)     2.5282   2.7194   2.9144
     
 M, MOL WT         24.147   23.708   23.235
 (DLV/DLP)T      -1.01404 -1.02343 -1.03346
 (DLV/DLT)P        1.2733   1.4747   1.7160
 CP, CAL/(G)(K)    0.8992   1.2693   1.7757
 GAMMA (S)         1.1552   1.1366   1.1215
 SON VEL,M/SEC     1157.2   1121.1   1082.8
     
 TRANSPORT PROPERTIES (GASES ONLY)
   CONDUCTIVITY IN UNITS OF MILLICALORIES/(CM)(K)(SEC)

 VISC,MILLIPOISE  0.95584  0.91493  0.86937
0  WITH EQUILIBRIUM REACTIONS

 CP, CAL/(G)(K)    0.8992   1.2693   1.7755
 CONDUCITVITY      1.6654   2.4794   3.5937
 PRANDTL NUMBER    0.5161   0.4684   0.4295
0  WITH FROZEN REACTIONS

 CP, CAL/(G)(K)    0.4517   0.4488   0.4456
 CONDUCITVITY      0.6947   0.6614   0.6262
 PRANDTL NUMBER    0.6215   0.6209   0.6187
     
0MOLE FRACTIONS

 CO               0.16668  0.17015  0.17172
 CO2              0.09703  0.08878  0.08205
 H                0.01116  0.02102  0.03429
 HCO RAD          0.00001  0.00000  0.00000
 HNO              0.00001  0.00000  0.00000
 HO2              0.00003  0.00002  0.00001
 H2               0.06776  0.07333  0.07838
 H2O              0.31097  0.28874  0.26649
 H2O2             0.00001  0.00000  0.00000
 N                0.00001  0.00001  0.00000
 NH               0.00001  0.00000  0.00000
 NH3              0.00001  0.00000  0.00000
 NO               0.00609  0.00711  0.00679
 N2               0.31168  0.30545  0.29946
 O                0.00216  0.00538  0.00994
 OH               0.02246  0.03157  0.03727
 O2               0.00390  0.00842  0.01358
0ADDITIONAL PRODUCTS WHICH WERE CONSIDERED BUT WHOSE MOLE FRACTIONS WERE LESS THAN  0.50000E-05 FOR ALL ASSIGNED CONDITIONS

  CH3NO             CH3ONO            HCO2              CH2NO2            CH3NO2            C                 CH              
  CH2               FORMALDEHYDE      FORMIC ACID       CH3               HYDROXYMETHYLENE  METHYLOXIDE       CH4             
  METHANOL          CN                NCN RAD           CNN RAD           C2                C2H RAD           ACETYLENE       
  KETENE            C2H3 RAD          METHYL CYANIDE    CH3CO RAD         CH2CHO RAD        ETHYLENE          ACETALDEHYDE    
  ACETIC ACID       (FORMIC ACID)2    ETHYL RAD         ETHYL OXIDE RAD   ETHANE            AZOMETHANE        DIMETHYL ETHER  
  ETHANOL           CNC RAD           CYANOGEN          CCO RAD           C3                C3H3 RAD          CYCLOPROPENE    
  PROPYNE           ALLENE            C3H5 RAD          CYCLOPROPANE      PROPYLENE         PROPYLENE OXIDE   I-PROPYL RAD    
  N-PROPYL RAD      PROPANE           1-PROPANOL        CARBON SUBOXIDE   C4                BUTADIYNE         CYCLOBUTADIENE  
  BUTAN-1EN-3YN     1,3-BUTADIENE     2-BUTYNE          2-BUTENE TRANS    2-BUTENE CIS      ISOBUTENE         1-BUTENE        
  (ACETIC ACID)2    T-BUTYL RAD       S-BUTYL RAD       N-BUTYL RAD       N-BUTANE          ISOBUTANE         CARBON SUBNITRID
  C5                CYCLOPENTADIENE   CYCLOPENTANE      1-PENTENE         T-PENTYL RAD      N-PENTYL RAD      PENTANE         
  ISOPENTANE        CH3C(CH3)2CH3     HEXATRIYNE        PHENYL RAD        PHENOXY RAD       BENZENE           PHENOL          
  CYCLOHEXENE       N-HEXYL RAD       TOLUENE           CRESOL            1-HEPTENE         N-HEPTYL RAD      N-HEPTANE       
  1-OCTENE          N-OCTYL RAD       OCTANE            ISO-OCTANE        N-NONYL RAD       NAPTHLENE         AZULENE         
  N-DECYL RAD       O-BIPHENYL RAD    BIPHENYL          JET-A(G)          HCN               HNCO              HNO2            
  HNO3              H2N2              NCO               NH2               NH2OH             NO2               NO3             
  N2H2              NH2NO2            N2H4              N2O               N2O3              N2O4              N2O5            
  N3                N3H               O3                C(GR)             METHANOL(L)       BENZENE(L)        TOLUENE(L)      
  OCTANE(L)         JET-A(L)          H2O(S)            H2O(L)          
0NOTE. WEIGHT FRACTION OF FUEL IN TOTAL FUELS AND OF OXIDANT IN TOTAL OXIDANTS
==================================
END: Output File: "cet86.out"
==================================

