

 ############
 # Title: [CurveBeam1]
 ############
 reset
 set terminal png small size  800, 800;
 set output 'FIGURES\CurveBeam1Qcv1Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "CurveBeam1";
 set xlabel "X1 [m]";
 set xrange [-6.00000e-002 :  6.00000e-002];
 set ylabel "X2 [m]";
 set yrange [ 0.00000e+000 :  1.20000e-001];
 plot 'FIGURES\CurveBeam1.mdt' using 1:2 title "Curve" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [CurveBeam1]
 ############
 reset
 set terminal png small size  800, 800;
 set output 'FIGURES\CurveBeam1Qcv2Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "CurveBeam1";
 set xlabel "X1 [m]";
 set xrange [-6.00000e-002 :  6.00000e-002];
 set ylabel "X3 [m]";
 set yrange [-6.00000e-002 :  6.00000e-002];
 plot 'FIGURES\CurveBeam1.mdt' using 1:3 title "Curve" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [CurveBeam1]
 ############
 reset
 set terminal png small size  800, 800;
 set output 'FIGURES\CurveBeam1Qcv3Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "CurveBeam1";
 set xlabel "X2 [m]";
 set xrange [ 0.00000e+000 :  1.20000e-001];
 set ylabel "X3 [m]";
 set yrange [-6.00000e-002 :  6.00000e-002];
 plot 'FIGURES\CurveBeam1.mdt' using 2:3 title "Curve" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [CurveBeam2]
 ############
 reset
 set terminal png small size  800, 800;
 set output 'FIGURES\CurveBeam2Qcv1Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "CurveBeam2";
 set xlabel "X1 [m]";
 set xrange [ 0.00000e+000 :  2.40000e-001];
 set ylabel "X2 [m]";
 set yrange [ 0.00000e+000 :  2.40000e-001];
 plot 'FIGURES\CurveBeam2.mdt' using 1:2 title "Curve" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [CurveBeam2]
 ############
 reset
 set terminal png small size  800, 800;
 set output 'FIGURES\CurveBeam2Qcv2Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "CurveBeam2";
 set xlabel "X1 [m]";
 set xrange [ 0.00000e+000 :  2.40000e-001];
 set ylabel "X3 [m]";
 set yrange [-1.20000e-001 :  1.20000e-001];
 plot 'FIGURES\CurveBeam2.mdt' using 1:3 title "Curve" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [CurveBeam2]
 ############
 reset
 set terminal png small size  800, 800;
 set output 'FIGURES\CurveBeam2Qcv3Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "CurveBeam2";
 set xlabel "X2 [m]";
 set xrange [ 0.00000e+000 :  2.40000e-001];
 set ylabel "X3 [m]";
 set yrange [-1.20000e-001 :  1.20000e-001];
 plot 'FIGURES\CurveBeam2.mdt' using 2:3 title "Curve" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [CurveBeam3]
 ############
 reset
 set terminal png small size  800, 800;
 set output 'FIGURES\CurveBeam3Qcv1Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "CurveBeam3";
 set xlabel "X1 [m]";
 set xrange [ 1.80000e-001 :  3.00000e-001];
 set ylabel "X2 [m]";
 set yrange [ 0.00000e+000 :  1.20000e-001];
 plot 'FIGURES\CurveBeam3.mdt' using 1:2 title "Curve" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [CurveBeam3]
 ############
 reset
 set terminal png small size  800, 800;
 set output 'FIGURES\CurveBeam3Qcv2Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "CurveBeam3";
 set xlabel "X1 [m]";
 set xrange [ 1.80000e-001 :  3.00000e-001];
 set ylabel "X3 [m]";
 set yrange [-6.00000e-002 :  6.00000e-002];
 plot 'FIGURES\CurveBeam3.mdt' using 1:3 title "Curve" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [CurveBeam3]
 ############
 reset
 set terminal png small size  800, 800;
 set output 'FIGURES\CurveBeam3Qcv3Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "CurveBeam3";
 set xlabel "X2 [m]";
 set xrange [ 0.00000e+000 :  1.20000e-001];
 set ylabel "X3 [m]";
 set yrange [-6.00000e-002 :  6.00000e-002];
 plot 'FIGURES\CurveBeam3.mdt' using 2:3 title "Curve" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [CurveTest]
 ############
 reset
 set terminal png small size  800, 800;
 set output 'FIGURES\CurveTestQcv1Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "CurveTest";
 set xlabel "X1 [m]";
 set xrange [ 0.00000e+000 :  1.00000e+000];
 set ylabel "X2 [m]";
 set yrange [-2.66306e-001 :  7.33694e-001];
 plot 'FIGURES\CurveTest.mdt' using 1:2 title "Curve" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [CurveTest]
 ############
 reset
 set terminal png small size  800, 800;
 set output 'FIGURES\CurveTestQcv2Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "CurveTest";
 set xlabel "X1 [m]";
 set xrange [ 0.00000e+000 :  1.00000e+000];
 set ylabel "X3 [m]";
 set yrange [-5.00000e-001 :  5.00000e-001];
 plot 'FIGURES\CurveTest.mdt' using 1:3 title "Curve" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [CurveTest]
 ############
 reset
 set terminal png small size  800, 800;
 set output 'FIGURES\CurveTestQcv3Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "CurveTest";
 set xlabel "X2 [m]";
 set xrange [-2.66306e-001 :  7.33694e-001];
 set ylabel "X3 [m]";
 set yrange [-5.00000e-001 :  5.00000e-001];
 plot 'FIGURES\CurveTest.mdt' using 2:3 title "Curve" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [OriDistBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\OriDistBeam1.png';
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title "OriDistBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "CONFORMAL ROTATIONS";
 set autoscale y;
 plot 'FIGURES\OriDistBeam1.mdt' using 1:2 title "C1" with lines linestyle 1 , \
      'FIGURES\OriDistBeam1.mdt' using 1:3 title "C2" with lines linestyle 2 , \
      'FIGURES\OriDistBeam1.mdt' using 1:4 title "C3" with lines linestyle 3 ;
 set nomultiplot;
 set output;

 ############
 # Title: [OriDistBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\OriDistBeam2.png';
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title "OriDistBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "CONFORMAL ROTATIONS";
 set autoscale y;
 plot 'FIGURES\OriDistBeam2.mdt' using 1:2 title "C1" with lines linestyle 1 , \
      'FIGURES\OriDistBeam2.mdt' using 1:3 title "C2" with lines linestyle 2 , \
      'FIGURES\OriDistBeam2.mdt' using 1:4 title "C3" with lines linestyle 3 ;
 set nomultiplot;
 set output;

 ############
 # Title: [OriDistBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\OriDistBeam3.png';
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title "OriDistBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "CONFORMAL ROTATIONS";
 set autoscale y;
 plot 'FIGURES\OriDistBeam3.mdt' using 1:2 title "C1" with lines linestyle 1 , \
      'FIGURES\OriDistBeam3.mdt' using 1:3 title "C2" with lines linestyle 2 , \
      'FIGURES\OriDistBeam3.mdt' using 1:4 title "C3" with lines linestyle 3 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1QaxsQ.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "AXIAL STIFFNESS [N]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1QaxsQ.mdt' using 1:2 title "axs" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1Qi22Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "BENDING STIFFNESS I22 [N.m^2]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1Qi22Q.mdt' using 1:2 title "i22" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1Qi33Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "BENDING STIFFNESS I33 [N.m^2]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1Qi33Q.mdt' using 1:2 title "i33" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1Qi23Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "BENDING STIFFNESS I23 [N.m^2]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1Qi23Q.mdt' using 1:2 title "i23" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1QtosQ.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "TORSIONAL STIFFNESS [N.m^2]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1QtosQ.mdt' using 1:2 title "tos" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1Qk22Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "SHEARING STIFFNESS K22 [N]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1Qk22Q.mdt' using 1:2 title "k22" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1Qk33Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "SHEARING STIFFNESS K33 [N]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1Qk33Q.mdt' using 1:2 title "k33" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1Qk23Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "SHEARING STIFFNESS K23 [N]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1Qk23Q.mdt' using 1:2 title "k23" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1Qm00Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "MASS/SPAN M00 [kg/m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1Qm00Q.mdt' using 1:2 title "m00" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1Qm11Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "MOMENT OF INERTIA/SPAN M11 [kg.m^2/m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1Qm11Q.mdt' using 1:2 title "m11" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1Qm22Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "MOMENT OF INERTIA/SPAN M22 [kg.m^2/m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1Qm22Q.mdt' using 1:2 title "m22" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1Qm33Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "MOMENT OF INERTIA/SPAN M33 [kg.m^2/m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1Qm33Q.mdt' using 1:2 title "m33" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1Qxm2Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "MASS CENTER LOCATION XM2 [m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1Qxm2Q.mdt' using 1:2 title "xm2" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1Qxm3Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "MASS CENTER LOCATION XM3 [m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1Qxm3Q.mdt' using 1:2 title "xm3" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1Qxk2Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "SHEAR CENTER LOCATION XK2 [m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1Qxk2Q.mdt' using 1:2 title "xk2" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1Qxk3Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "SHEAR CENTER LOCATION XK3 [m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1Qxk3Q.mdt' using 1:2 title "xk3" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1Qxc2Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "CENTROID LOCATION XC2 [m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1Qxc2Q.mdt' using 1:2 title "xc2" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1Qxc3Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "CENTROID LOCATION XC3 [m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1Qxc3Q.mdt' using 1:2 title "xc3" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam1QmucQ.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam1";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "DAMPING COEFFICIENT [1/sec]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam1QmucQ.mdt' using 1:2 title "muc" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2QaxsQ.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "AXIAL STIFFNESS [N]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2QaxsQ.mdt' using 1:2 title "axs" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2Qi22Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "BENDING STIFFNESS I22 [N.m^2]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2Qi22Q.mdt' using 1:2 title "i22" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2Qi33Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "BENDING STIFFNESS I33 [N.m^2]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2Qi33Q.mdt' using 1:2 title "i33" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2Qi23Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "BENDING STIFFNESS I23 [N.m^2]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2Qi23Q.mdt' using 1:2 title "i23" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2QtosQ.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "TORSIONAL STIFFNESS [N.m^2]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2QtosQ.mdt' using 1:2 title "tos" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2Qk22Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "SHEARING STIFFNESS K22 [N]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2Qk22Q.mdt' using 1:2 title "k22" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2Qk33Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "SHEARING STIFFNESS K33 [N]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2Qk33Q.mdt' using 1:2 title "k33" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2Qk23Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "SHEARING STIFFNESS K23 [N]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2Qk23Q.mdt' using 1:2 title "k23" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2Qm00Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "MASS/SPAN M00 [kg/m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2Qm00Q.mdt' using 1:2 title "m00" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2Qm11Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "MOMENT OF INERTIA/SPAN M11 [kg.m^2/m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2Qm11Q.mdt' using 1:2 title "m11" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2Qm22Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "MOMENT OF INERTIA/SPAN M22 [kg.m^2/m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2Qm22Q.mdt' using 1:2 title "m22" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2Qm33Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "MOMENT OF INERTIA/SPAN M33 [kg.m^2/m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2Qm33Q.mdt' using 1:2 title "m33" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2Qxm2Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "MASS CENTER LOCATION XM2 [m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2Qxm2Q.mdt' using 1:2 title "xm2" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2Qxm3Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "MASS CENTER LOCATION XM3 [m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2Qxm3Q.mdt' using 1:2 title "xm3" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2Qxk2Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "SHEAR CENTER LOCATION XK2 [m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2Qxk2Q.mdt' using 1:2 title "xk2" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2Qxk3Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "SHEAR CENTER LOCATION XK3 [m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2Qxk3Q.mdt' using 1:2 title "xk3" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2Qxc2Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "CENTROID LOCATION XC2 [m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2Qxc2Q.mdt' using 1:2 title "xc2" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2Qxc3Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "CENTROID LOCATION XC3 [m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2Qxc3Q.mdt' using 1:2 title "xc3" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam2]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam2QmucQ.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam2";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "DAMPING COEFFICIENT [1/sec]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam2QmucQ.mdt' using 1:2 title "muc" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3QaxsQ.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "AXIAL STIFFNESS [N]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3QaxsQ.mdt' using 1:2 title "axs" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3Qi22Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "BENDING STIFFNESS I22 [N.m^2]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3Qi22Q.mdt' using 1:2 title "i22" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3Qi33Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "BENDING STIFFNESS I33 [N.m^2]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3Qi33Q.mdt' using 1:2 title "i33" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3Qi23Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "BENDING STIFFNESS I23 [N.m^2]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3Qi23Q.mdt' using 1:2 title "i23" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3QtosQ.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "TORSIONAL STIFFNESS [N.m^2]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3QtosQ.mdt' using 1:2 title "tos" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3Qk22Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "SHEARING STIFFNESS K22 [N]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3Qk22Q.mdt' using 1:2 title "k22" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3Qk33Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "SHEARING STIFFNESS K33 [N]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3Qk33Q.mdt' using 1:2 title "k33" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3Qk23Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "SHEARING STIFFNESS K23 [N]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3Qk23Q.mdt' using 1:2 title "k23" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3Qm00Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "MASS/SPAN M00 [kg/m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3Qm00Q.mdt' using 1:2 title "m00" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3Qm11Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "MOMENT OF INERTIA/SPAN M11 [kg.m^2/m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3Qm11Q.mdt' using 1:2 title "m11" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3Qm22Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "MOMENT OF INERTIA/SPAN M22 [kg.m^2/m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3Qm22Q.mdt' using 1:2 title "m22" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3Qm33Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "MOMENT OF INERTIA/SPAN M33 [kg.m^2/m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3Qm33Q.mdt' using 1:2 title "m33" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3Qxm2Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "MASS CENTER LOCATION XM2 [m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3Qxm2Q.mdt' using 1:2 title "xm2" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3Qxm3Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "MASS CENTER LOCATION XM3 [m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3Qxm3Q.mdt' using 1:2 title "xm3" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3Qxk2Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "SHEAR CENTER LOCATION XK2 [m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3Qxk2Q.mdt' using 1:2 title "xk2" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3Qxk3Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "SHEAR CENTER LOCATION XK3 [m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3Qxk3Q.mdt' using 1:2 title "xk3" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3Qxc2Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "CENTROID LOCATION XC2 [m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3Qxc2Q.mdt' using 1:2 title "xc2" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3Qxc3Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "CENTROID LOCATION XC3 [m]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3Qxc3Q.mdt' using 1:2 title "xc3" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam3]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\PropertyBeam3QmucQ.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "PropertyBeam3";
 set xlabel "ETA";
 set autoscale x;
 set ylabel "DAMPING COEFFICIENT [1/sec]";
 set autoscale y;
 plot 'FIGURES\PropertyBeam3QmucQ.mdt' using 1:2 title "muc" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [ScheduleRotationA]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\ScheduleRotationA.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "ScheduleRotationA";
 set xlabel "TIME [sec]";
 set autoscale x;
 set ylabel "TIME FUNCTION";
 set autoscale y;
 plot 'FIGURES\ScheduleRotationA.mdt' using 1:2 title "History" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [ScheduleA]
 ############
 reset
 set terminal png small size  800, 600;
 set output 'FIGURES\ScheduleA.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "ScheduleA";
 set xlabel "TIME [sec]";
 set autoscale x;
 set ylabel "TIME FUNCTION";
 set autoscale y;
 plot 'FIGURES\ScheduleA.mdt' using 1:2 title "History" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [CurveEigthCircle]
 ############
 reset
 set terminal png small size  800, 800;
 set output 'FIGURES\CurveEigthCircleQcv1Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "CurveEigthCircle";
 set xlabel "X1 [m]";
 set xrange [-3.53553e-001 :  3.53553e-001];
 set ylabel "X2 [m]";
 set yrange [ 5.00000e-001 :  1.20711e+000];
 plot 'FIGURES\CurveEigthCircle.mdt' using 1:2 title "Curve" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [CurveEigthCircle]
 ############
 reset
 set terminal png small size  800, 800;
 set output 'FIGURES\CurveEigthCircleQcv2Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "CurveEigthCircle";
 set xlabel "X1 [m]";
 set xrange [-3.53553e-001 :  3.53553e-001];
 set ylabel "X3 [m]";
 set yrange [ 0.00000e+000 :  7.07107e-001];
 plot 'FIGURES\CurveEigthCircle.mdt' using 1:3 title "Curve" with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [CurveEigthCircle]
 ############
 reset
 set terminal png small size  800, 800;
 set output 'FIGURES\CurveEigthCircleQcv3Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title "CurveEigthCircle";
 set xlabel "X2 [m]";
 set xrange [ 5.00000e-001 :  1.20711e+000];
 set ylabel "X3 [m]";
 set yrange [ 0.00000e+000 :  7.07107e-001];
 plot 'FIGURES\CurveEigthCircle.mdt' using 2:3 title "Curve" with lines linestyle 1 ;
 set nomultiplot;
 set output;