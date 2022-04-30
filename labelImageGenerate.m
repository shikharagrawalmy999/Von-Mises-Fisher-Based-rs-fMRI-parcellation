function [outputLabelImage, outputTimeSeries, outputNumberOfLabels] = labelImageGenerate(n, l, rho_in, rho_ac, T)
    N = n;
    upperLimit = 20; % As was mentioned in the srikanth et. al
    L = l;
    beta_s = 80;
    %beta_s = 2;
    image = upperLimit*rand([N N]);
    means = (0:L-1);
    means = means/(L-1);
    means = means + 1/(2*(L-1));
    means = upperLimit*means;
    means = means';
    stdDev = ones(L,1);
    var = stdDev.*stdDev;
    labelImage = randi([1, L],N,N);
    count = zeros(L,1);
    for m=1:50
        for i=1:N
            for j=1:N
                for l=1:L
                    count(l) = (labelImage(mod(i+1, N) + N*(mod(i+1, N)==0), mod(j, N) + N*(mod(j, N)==0)) ~= l)+(labelImage(mod(i-1, N) + N*(mod(i-1, N)==0), mod(j, N) + N*(mod(j, N)==0)) ~= l)+(labelImage(mod(i, N) + N*(mod(i, N)==0), mod(j+1, N) + N*(mod(j+1, N)==0)) ~= l)+(labelImage(mod(i, N) + N*(mod(i, N)==0), mod(j-1, N) + N*(mod(j-1, N)==0)) ~= l); 
                end
                value = (exp(-((image(i,j)-means).^2)./(2*var)).*exp(-(count*beta_s)))./(2*pi*stdDev);
                [M,I]=max(value);
                labelImage(i,j)=I;
            end
        end
    end
    outputLabelImage = labelImage;
    imagesc(labelImage);
    cov = zeros(N*N, N*N);
    for i=1:N*N
        for j=1:N*N
            if i==j
                cov(i,j) = 1;
            else
                row1 = floor((i-1)/N) + 1;
                column1 = mod(i-1,N)+1;
                row2 = floor((j-1)/N) + 1;
                column2 = mod(j-1,N)+1;
                cov(i,j) = rho_in*(labelImage(row1, column1)==labelImage(row2,column2))+rho_ac*(labelImage(row1, column1)~=labelImage(row2,column2));
            end
        end
    end
    [V,D] = eig(cov);
    A = V*sqrt(D);
    timeSeries = zeros(T, N*N);
    for i=1:T
        timeSeries(i, :) = (A*randn(N*N, 1))';
    end
    norms = sqrt(sum(timeSeries.*timeSeries));
    norms = repmat(norms, T, 1);
    timeSeries = timeSeries./norms;
    outputTimeSeries = timeSeries;
    
    temp = size(unique(reshape(labelImage, [N*N, 1])));
    outputNumberOfLabels = temp(1);
    
end