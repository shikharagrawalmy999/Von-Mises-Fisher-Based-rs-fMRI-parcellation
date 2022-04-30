L = 5;
N = 60;
T = 20;
rho_in = 0.9;
rho_ac = 0.1;
[realLabelImage, timeSeries, realNumberOfLabels] = labelImageGenerate(N, L, rho_in, rho_ac, T);