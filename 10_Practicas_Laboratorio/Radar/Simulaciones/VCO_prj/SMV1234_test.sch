<Qucs Schematic 0.0.20>
<Properties>
  <View=50,-7,1182,600,0.90636,0,0>
  <Grid=10,10,1>
  <DataSet=SMV1234_test.dat>
  <DataDisplay=SMV1234_test.dpl>
  <OpenDisplay=1>
  <Script=SMV1234_test.m>
  <RunScript=0>
  <showFrame=0>
  <FrameText0=Título>
  <FrameText1=Dibujado por:>
  <FrameText2=Fecha:>
  <FrameText3=Revisión:>
</Properties>
<Symbol>
</Symbol>
<Components>
  <.DC DC1 1 100 310 0 39 0 0 "26.85" 0 "0.001" 0 "1 pA" 0 "1 uV" 0 "no" 0 "150" 0 "no" 0 "none" 0 "CroutLU" 0>
  <Pac P2 1 380 420 18 -26 0 1 "1" 1 "50 Ohm" 1 "0 dBm" 0 "1 GHz" 0 "26.85" 0>
  <GND *5 5 380 490 0 0 0 0>
  <.SP SP1 1 100 400 0 63 0 0 "lin" 1 "100 MHz" 1 "1 GHz" 1 "1000" 1 "no" 0 "1" 0 "2" 0 "no" 0 "no" 0>
  <GND *1 5 690 430 0 0 0 0>
  <Pac P3 1 820 420 18 -26 0 1 "2" 1 "50 Ohm" 1 "0 dBm" 0 "1 GHz" 0 "26.85" 0>
  <GND *6 5 820 490 0 0 0 0>
  <GND *7 5 1130 430 0 0 0 0>
  <GND * 5 1050 390 0 0 0 0>
  <SPfile X1 1 1050 360 -26 -59 0 0 "/home/aalmela/scm/git/utn/electronicaaplicada3/10_Practicas_Laboratorio/RadarDoppler900MHZ/LNA/simulaciones/QUCS/LNA_radar_prj/capacitor/sparameter-mlcc20221002200639/grm-s-v60/Temperature_compensating/2012_0805/GRM21A5C2E470FW01.s2p" 0 "rectangular" 0 "linear" 0 "open" 0 "2" 0>
  <Sub SUB1 1 610 360 -26 21 1 2 "GRM21A5C2E470FW01_spice.sch" 0>
</Components>
<Wires>
  <640 360 690 360 "" 0 0 0 "">
  <380 450 380 490 "" 0 0 0 "">
  <690 360 690 430 "" 0 0 0 "">
  <380 360 380 390 "" 0 0 0 "">
  <380 360 580 360 "" 0 0 0 "">
  <1080 360 1130 360 "" 0 0 0 "">
  <820 450 820 490 "" 0 0 0 "">
  <1130 360 1130 430 "" 0 0 0 "">
  <820 360 820 390 "" 0 0 0 "">
  <820 360 1020 360 "" 0 0 0 "">
</Wires>
<Diagrams>
  <Smith 350 290 257 257 3 #c0c0c0 1 00 1 0 1 1 1 0 4 1 1 0 1 1 315 0 225 "" "" "" "">
	<"S[1,1]" #0000ff 0 3 0 0 0>
	  <Mkr 9.66667e+08/0/0 -239 -243 3 0 0 10>
	  <Mkr 8e+08/0/0 -239 -85 3 0 0 10>
	  <Mkr 1e+08/0/0 53 -203 3 0 0 10>
  </Smith>
  <Smith 870 280 257 257 3 #c0c0c0 1 00 1 0 1 1 1 0 4 1 1 0 1 1 315 0 225 "" "" "" "">
	<"S[2,2]" #ff0000 0 3 0 0 0>
	  <Mkr 9.66667e+08/0/0 -199 -253 3 0 0 10>
	  <Mkr 8e+08/0/0 -189 -75 3 0 0 10>
	  <Mkr 1.03604e+08/0/0 56 -152 3 0 0 10>
  </Smith>
</Diagrams>
<Paintings>
</Paintings>
