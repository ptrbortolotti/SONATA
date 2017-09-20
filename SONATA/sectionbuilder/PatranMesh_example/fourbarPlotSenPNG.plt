

 ############
 # Title: [SensorBeam1Forces]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\SensorBeam1ForcesQforQ.png';
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title "SensorBeam1Forces";
 set xlabel "TIME [sec]";
 set autoscale x;
 set ylabel "FORCE [N]";
 set autoscale y;
 plot 'FIGURES\SensorBeam1Forces.mdt' using 1:2 title "F_1" with lines linestyle 1 , \
      'FIGURES\SensorBeam1Forces.mdt' using 1:3 title "F_2" with lines linestyle 2 , \
      'FIGURES\SensorBeam1Forces.mdt' using 1:4 title "F_3" with lines linestyle 3 ;
 set nomultiplot;
 set output;

 ############
 # Title: [SensorBeam1Forces]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\SensorBeam1ForcesQmomQ.png';
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title "SensorBeam1Forces";
 set xlabel "TIME [sec]";
 set autoscale x;
 set ylabel "MOMENT [N.m]";
 set autoscale y;
 plot 'FIGURES\SensorBeam1Forces.mdt' using 1:5 title "M_1" with lines linestyle 1 , \
      'FIGURES\SensorBeam1Forces.mdt' using 1:6 title "M_2" with lines linestyle 2 , \
      'FIGURES\SensorBeam1Forces.mdt' using 1:7 title "M_3" with lines linestyle 3 ;
 set nomultiplot;
 set output;

 ############
 # Title: [SensorBeam1Velocities]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\SensorBeam1VelocitiesQvelQ.png';
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title "SensorBeam1Velocities";
 set xlabel "TIME [sec]";
 set autoscale x;
 set ylabel "VELOCITY [m/sec]";
 set autoscale y;
 plot 'FIGURES\SensorBeam1Velocities.mdt' using 1:2 title "V_1" with lines linestyle 1 , \
      'FIGURES\SensorBeam1Velocities.mdt' using 1:3 title "V_2" with lines linestyle 2 , \
      'FIGURES\SensorBeam1Velocities.mdt' using 1:4 title "V_3" with lines linestyle 3 ;
 set nomultiplot;
 set output;

 ############
 # Title: [SensorBeam1Velocities]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\SensorBeam1VelocitiesQomeQ.png';
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title "SensorBeam1Velocities";
 set xlabel "TIME [sec]";
 set autoscale x;
 set ylabel "ANGULAR VELOCITY [rad/sec]";
 set autoscale y;
 plot 'FIGURES\SensorBeam1Velocities.mdt' using 1:5 title "W_1" with lines linestyle 1 , \
      'FIGURES\SensorBeam1Velocities.mdt' using 1:6 title "W_2" with lines linestyle 2 , \
      'FIGURES\SensorBeam1Velocities.mdt' using 1:7 title "W_3" with lines linestyle 3 ;
 set nomultiplot;
 set output;

 ############
 # Title: [SensorBeam2Displacements]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\SensorBeam2DisplacementsQdisQ.png';
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title "SensorBeam2Displacements";
 set xlabel "TIME [sec]";
 set autoscale x;
 set ylabel "DISPLACEMENTS [m]";
 set autoscale y;
 plot 'FIGURES\SensorBeam2Displacements.mdt' using 1:2 title "u_1" with lines linestyle 1 , \
      'FIGURES\SensorBeam2Displacements.mdt' using 1:3 title "u_2" with lines linestyle 2 , \
      'FIGURES\SensorBeam2Displacements.mdt' using 1:4 title "u_3" with lines linestyle 3 ;
 set nomultiplot;
 set output;

 ############
 # Title: [SensorBeam2Displacements]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\SensorBeam2DisplacementsQrotQ.png';
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title "SensorBeam2Displacements";
 set xlabel "TIME [sec]";
 set autoscale x;
 set ylabel "EULER ANGLES (3-1-2 SEQUENCE) [deg]";
 set autoscale y;
 plot 'FIGURES\SensorBeam2Displacements.mdt' using 1:5 title "Phi" with lines linestyle 1 , \
      'FIGURES\SensorBeam2Displacements.mdt' using 1:6 title "Teta" with lines linestyle 2 , \
      'FIGURES\SensorBeam2Displacements.mdt' using 1:7 title "Psi" with lines linestyle 3 ;
 set nomultiplot;
 set output;

 ############
 # Title: [SensorBeam2Forces]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\SensorBeam2ForcesQforQ.png';
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title "SensorBeam2Forces";
 set xlabel "TIME [sec]";
 set autoscale x;
 set ylabel "FORCE [N]";
 set autoscale y;
 plot 'FIGURES\SensorBeam2Forces.mdt' using 1:2 title "F_1" with lines linestyle 1 , \
      'FIGURES\SensorBeam2Forces.mdt' using 1:3 title "F_2" with lines linestyle 2 , \
      'FIGURES\SensorBeam2Forces.mdt' using 1:4 title "F_3" with lines linestyle 3 ;
 set nomultiplot;
 set output;

 ############
 # Title: [SensorBeam2Forces]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\SensorBeam2ForcesQmomQ.png';
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title "SensorBeam2Forces";
 set xlabel "TIME [sec]";
 set autoscale x;
 set ylabel "MOMENT [N.m]";
 set autoscale y;
 plot 'FIGURES\SensorBeam2Forces.mdt' using 1:5 title "M_1" with lines linestyle 1 , \
      'FIGURES\SensorBeam2Forces.mdt' using 1:6 title "M_2" with lines linestyle 2 , \
      'FIGURES\SensorBeam2Forces.mdt' using 1:7 title "M_3" with lines linestyle 3 ;
 set nomultiplot;
 set output;

 ############
 # Title: [SensorBeam2Velocities]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\SensorBeam2VelocitiesQvelQ.png';
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title "SensorBeam2Velocities";
 set xlabel "TIME [sec]";
 set autoscale x;
 set ylabel "VELOCITY [m/sec]";
 set autoscale y;
 plot 'FIGURES\SensorBeam2Velocities.mdt' using 1:2 title "V_1" with lines linestyle 1 , \
      'FIGURES\SensorBeam2Velocities.mdt' using 1:3 title "V_2" with lines linestyle 2 , \
      'FIGURES\SensorBeam2Velocities.mdt' using 1:4 title "V_3" with lines linestyle 3 ;
 set nomultiplot;
 set output;

 ############
 # Title: [SensorBeam2Velocities]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\SensorBeam2VelocitiesQomeQ.png';
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title "SensorBeam2Velocities";
 set xlabel "TIME [sec]";
 set autoscale x;
 set ylabel "ANGULAR VELOCITY [rad/sec]";
 set autoscale y;
 plot 'FIGURES\SensorBeam2Velocities.mdt' using 1:5 title "W_1" with lines linestyle 1 , \
      'FIGURES\SensorBeam2Velocities.mdt' using 1:6 title "W_2" with lines linestyle 2 , \
      'FIGURES\SensorBeam2Velocities.mdt' using 1:7 title "W_3" with lines linestyle 3 ;
 set nomultiplot;
 set output;

 ############
 # Title: [SensorRotationA]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\SensorRotationA.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "SensorRotationA";
 set xlabel "TIME [sec]";
 set autoscale x;
 set ylabel "RELATIVE ROTATION [rad]";
 set autoscale y;
 plot 'FIGURES\SensorRotationA.mdt' using 1:5 title "phi" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [SensorRotationB]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\SensorRotationB.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "SensorRotationB";
 set xlabel "TIME [sec]";
 set autoscale x;
 set ylabel "RELATIVE ROTATION [rad]";
 set autoscale y;
 plot 'FIGURES\SensorRotationB.mdt' using 1:5 title "phi" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [SensorRotationC]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\SensorRotationC.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "SensorRotationC";
 set xlabel "TIME [sec]";
 set autoscale x;
 set ylabel "RELATIVE ROTATION [rad]";
 set autoscale y;
 plot 'FIGURES\SensorRotationC.mdt' using 1:5 title "phi" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [SensorRotationD]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\SensorRotationD.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "SensorRotationD";
 set xlabel "TIME [sec]";
 set autoscale x;
 set ylabel "RELATIVE ROTATION [rad]";
 set autoscale y;
 plot 'FIGURES\SensorRotationD.mdt' using 1:5 title "phi" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [SensorPrdA]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\SensorPrdAQprdrQ.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "SensorPrdA";
 set xlabel "TIME [sec]";
 set autoscale x;
 set ylabel "PRESCRIBED ROTATION [rad]";
 set autoscale y;
 plot 'FIGURES\SensorPrdA.mdt' using 1:2 title "phi" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [SensorPrdA]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\SensorPrdAQprdmQ.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "SensorPrdA";
 set xlabel "TIME [sec]";
 set autoscale x;
 set ylabel "DRIVING MOMENT [N.m]";
 set autoscale y;
 plot 'FIGURES\SensorPrdA.mdt' using 1:3 title "Md" with lines linestyle 1 ;
 set nomultiplot;
 set output;