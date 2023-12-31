<Qucs Schematic 0.0.20>
<Properties>
  <View=306,119,1583,704,1.4852,0,19>
  <Grid=10,10,1>
  <DataSet=EsquemaLNA.dat>
  <DataDisplay=EsquemaLNA.dpl>
  <OpenDisplay=1>
  <Script=EsquemaLNA.m>
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
  <GND * 5 400 380 0 0 0 0>
  <GND * 5 740 330 0 0 0 0>
  <L L1 1 490 300 -26 -56 1 0 "100 nH" 1 "" 0>
  <R R1 1 700 220 -26 -61 1 0 "1800 Ohm" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "US" 0>
  <C C1 1 780 220 -26 17 0 0 "330 pF" 1 "" 0 "neutral" 0>
  <R R2 1 850 450 -93 -26 1 1 "68 Ohm" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "US" 0>
  <GND * 5 850 580 0 0 0 0>
  <GND * 5 1110 390 0 0 0 0>
  <L L2 1 1040 300 -26 10 0 0 "39 nH" 1 "" 0>
  <GND * 5 970 450 0 0 0 0>
  <C C6 1 910 300 -26 17 0 0 "330 pF" 1 "" 0 "neutral" 0>
  <R R3 1 650 400 15 -26 0 1 "47k Ohm" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "US" 0>
  <C C5 1 850 530 17 -26 0 1 "47 nF" 1 "" 0 "neutral" 0>
  <GND * 5 540 380 0 0 0 0>
  <C C7 1 590 300 -26 -63 1 0 "330p" 1 "" 0 "neutral" 0>
  <Eqn Eqn1 1 960 200 -32 19 0 0 "dBS21=dB(S[2,1])" 1 "yes" 0>
  <Pac P1 1 400 330 18 -26 0 1 "1" 1 "50 Ohm" 1 "-40 dBm" 0 "1 GHz" 0 "26.85" 0>
  <Pac P2 1 1110 340 18 -26 0 1 "2" 1 "50 Ohm" 1 "-40 dBm" 0 "1 GHz" 0 "26.85" 0>
  <SPfile X1 1 740 300 -26 -43 0 0 "/home/aalmela/scm/git/utn/eaiii/2019/proyecto/LNA/modelos/SPAR/BFP460/BFP460_w_noise_VCE_2.5V_IC_4.0mA.s2p" 0 "rectangular" 0 "linear" 0 "open" 0 "2" 0>
  <C C2 1 540 330 17 -26 1 3 "18 pF" 1 "" 0 "neutral" 0>
  <C C4 1 970 380 17 -26 0 1 "10 pF" 1 "" 0 "neutral" 0>
  <.SP SP1 1 460 470 0 79 0 0 "lin" 1 "70 MHz" 1 "130 MHz" 1 "101" 1 "yes" 0 "1" 0 "2" 0 "no" 0 "no" 0>
</Components>
<Wires>
  <400 360 400 380 "" 0 0 0 "">
  <400 300 460 300 "" 0 0 0 "">
  <650 300 710 300 "" 0 0 0 "">
  <650 220 650 300 "" 0 0 0 "">
  <650 220 670 220 "" 0 0 0 "">
  <730 220 750 220 "" 0 0 0 "">
  <770 300 850 300 "" 0 0 0 "">
  <520 300 540 300 "" 0 0 0 "">
  <850 480 850 490 "" 0 0 0 "">
  <850 560 850 580 "" 0 0 0 "">
  <850 300 850 420 "" 0 0 0 "">
  <1110 300 1110 310 "" 0 0 0 "">
  <1110 370 1110 390 "" 0 0 0 "">
  <1070 300 1110 300 "" 0 0 0 "">
  <850 300 880 300 "" 0 0 0 "">
  <970 300 1010 300 "" 0 0 0 "">
  <970 300 970 350 "" 0 0 0 "">
  <970 410 970 450 "" 0 0 0 "">
  <940 300 970 300 "" 0 0 0 "">
  <650 300 650 370 "" 0 0 0 "">
  <650 430 650 490 "" 0 0 0 "">
  <850 490 850 500 "" 0 0 0 "">
  <650 490 850 490 "" 0 0 0 "">
  <540 360 540 380 "" 0 0 0 "">
  <620 300 650 300 "" 0 0 0 "">
  <810 220 850 220 "" 0 0 0 "">
  <850 220 850 300 "" 0 0 0 "">
  <540 300 560 300 "" 0 0 0 "">
</Wires>
<Diagrams>
</Diagrams>
<Paintings>
</Paintings>
