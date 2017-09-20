

 ############
 # Title: [SurveyBeam2Displacements]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\SurveyBeam2DisplacementsQdisQ.png';
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title "SurveyBeam2Displacements";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "DISPLACEMENTS [m]";
 set autoscale y;
 plot 'FIGURES\SurveyBeam2Displacements.mdt' using 1:2 title "u_1" with lines linestyle 1 , \
      'FIGURES\SurveyBeam2Displacements.mdt' using 1:3 title "u_2" with lines linestyle 2 , \
      'FIGURES\SurveyBeam2Displacements.mdt' using 1:4 title "u_3" with lines linestyle 3 ;
 set nomultiplot;
 set output;

 ############
 # Title: [SurveyBeam2Displacements]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\SurveyBeam2DisplacementsQrotQ.png';
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title "SurveyBeam2Displacements";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "ROTATIONS";
 set autoscale y;
 plot 'FIGURES\SurveyBeam2Displacements.mdt' using 1:5 title "R_1" with lines linestyle 1 , \
      'FIGURES\SurveyBeam2Displacements.mdt' using 1:6 title "R_2" with lines linestyle 2 , \
      'FIGURES\SurveyBeam2Displacements.mdt' using 1:7 title "R_3" with lines linestyle 3 ;
 set nomultiplot;
 set output;