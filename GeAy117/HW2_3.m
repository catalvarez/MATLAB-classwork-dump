[ lowerbound0568,upperbound0568,peak0568 ] = HW2_1D_confidence(H05,P05,68);
[ lowerbound0595,upperbound0595,peak0595 ] = HW2_1D_confidence(H05,P05,95);
[ lowerbound1068,upperbound1068,peak1068 ] = HW2_1D_confidence(H10,P10,68);
[ lowerbound1095,upperbound1095,peak1095 ] = HW2_1D_confidence(H10,P10,95);
[ lowerbound2068,upperbound2068,peak2068 ] = HW2_1D_confidence(H20,P20,68);
[ lowerbound2095,upperbound2095,peak2095 ] = HW2_1D_confidence(H20,P20,95);
[ lowerbound5068,upperbound5068,peak5068 ] = HW2_1D_confidence(H50,P50,68);
[ lowerbound5095,upperbound5095,peak5095 ] = HW2_1D_confidence(H50,P50,95);

array = [lowerbound0568,upperbound0568,lowerbound0595,upperbound0595,peak0595;...
    lowerbound1068,upperbound1068,lowerbound1095,upperbound1095,peak1095;...
    lowerbound2068,upperbound2068,lowerbound2095,upperbound2095,peak2095;...
    lowerbound5068,upperbound5068,lowerbound5095,upperbound5095,peak5095];

N =[5;10;20;50];
T = 0.5.*(array(:,2)-array(:,1));
E = sqrt(array(:,5).*(1-array(:,5))./N);
percenterror = (E-T)./T.*100;
plot(N,percenterror)
title('Gaussian vs True Error')
xlabel('Nimber of Coin Flips')
ylabel('Percent Error')