% res = [0.1 0.01 0.001];
% DOF8 = [ 3.052467e-08 2.297e-5 0.0023];
% SNR9 = [1 3 10 30 100 1000 3000 10000];
% DOF9 = [2.2976e+03 28.366 0.2298 0.0028 2.297e-5 2.2977e-09 2.8367e-11 2.2984e-13];
res = [1 0.1 0.01 0.001] ;
for i = 1:4
     [H,DOF]= RadiationHW4_3(res(i),100);
     H8(i) = H;
     DOF8(i) = DOF;
end
SNR = [1 3 10 30 100 300 1000 3000 10000 30000];
for i = 1:10
     [H,DOF]= RadiationHW4_3(0.01,SNR(i));
     H9(i) = H;
     DOF9(i) = DOF;
end
figure (1)
semilogx(res,DOF8)
title('Problem 8')
xlabel('Resolution')
ylabel('DOF')
figure (2)
semilogx(res,H8)
title('Problem 8')
xlabel('Resolution')
ylabel('H')
figure (3)
semilogx(SNR,DOF9)
title('Problem 9')
xlabel('SNR')
ylabel('DOF')
figure (4)
semilogx(SNR,H9)
title('Problem 9')
xlabel('SNR')
ylabel('H')