%Code to conduct VMF-MRF based parcellation using alpha-expansion
%Needs as input resting state fMRI time series
beta_s = 0.1;
beta_l = 0.1;
N = 60;
T = 20;
%Step 0: calculate random cluster initialization
numOfLabels = 10;
initialLabels = randi([1, numOfLabels],N,N);
setOfLabels = unique(reshape(initialLabels, [N^2, 1]));
numOfLabels = size(setOfLabels);
numOfLabels = numOfLabels(1);
setOfLabels = setOfLabels';
%Step 1: calculate objective function
obj_func = 10000000;
old_obj_func = 11000000;
kappa = zeros(10, 1);
%iterative algorithm to find parameter and label updates
for x=1:1
    %Step 2.1: do closed form parameter updates
    mean = zeros(T, 20);
    for i=1:numOfLabels
        s = 0;
        thisLabel = setOfLabels(i);
        for j=1:N^2
            if initialLabels(floor((j-1)/N) + 1, mod(j-1,N)+1)==thisLabel
                mean(:, thisLabel) = mean(:, thisLabel) + timeSeries(:, j);
                s = s + 1;
            end
        end
        
        if s==0 
            mean(:, thisLabel) = zeros(T, 1);
        else
            mean(:, thisLabel) = mean(:, thisLabel)/s;
        end
        
        r = norm(mean(:, thisLabel));
        %disp("I am r");
        %disp(r);
        if r > 0
           mean(:, i) = mean(:, i)/r; 
        end
        if abs(r - 1) < 0.001 
            kappa(i) = 100; % setting kappa to a large value if r==1
        else
            kappa(i) = (r*T - r^3)/(1 - r*r);
        end
        
        if kappa(i) < 0
            kappa(i) = 0.001;
        end
        
        %disp(kappa(i));
    end
    %Step 2.2: do alpha-expansion on labels
    timeSeries = reshape(timeSeries, [T, N, N]);
    tempSet = setOfLabels;
    for i=1:numOfLabels
        alpha = tempSet(i);
        setOfLabels
        if (any(setOfLabels(:) == alpha))
            [initialLabels, setOfLabels] = AlphaExpansion(timeSeries, initialLabels, alpha, setOfLabels, beta_s, beta_l, N, mean, kappa, T);
        end
    end
    timeSeries = reshape(timeSeries, [T, N*N]);
    numOfLabels = size(setOfLabels);
    numOfLabels = numOfLabels(2);
    %Step 2.3: Upadte objective function
    old_obj_func = obj_func;
    obj_func = energy(N, initialLabels, beta_s, beta_l, mean, kappa, timeSeries);
end

%Step 3: Give clusters as output