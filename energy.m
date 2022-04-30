function [outputEnergy] = energy(N, labels, beta_s, beta_l, mean, kappa, timeSeries)
    sz = size(timeSeries);
    T = sz(1);
    C_kappa = kappa.^(T/2 - 1);
    C_kappa = C_kappa./(besseli(T/2 - 1,kappa));
    C_kappa = C_kappa/((2*pi)^(T/2));
    
    mat1 = circshift(labels,1,2);
    mat2 = circshift(labels,-1,2);
    mat3 = circshift(labels,1,1);
    mat4 = circshift(labels,-1,1);
    
    mat1 = mat1 ~= labels;
    mat2 = mat2 ~= labels;
    mat3 = mat3 ~= labels;
    mat4 = mat4 ~= labels;
    
    penaltyMatrix = mat1 + mat2 + mat3 + mat4;
    % divide by 2 since each clique is considered twice in the above
    % computation
    penalty = sum(sum(penaltyMatrix))/2;
    
    % computing prior 
    log_prior_spatial = (beta_s*penalty);
    
    temp = size(unique(reshape(labels, [N*N, 1])));
    log_prior_label = (temp(1)*temp(2))*beta_l;
    log_likelihood = 0;
    
    for i=1:N
        for j=1:N
            %log_likelihood = log_likelihood + log(VMFMeanDirDensity((mean(:, labels(i,j)))'*timeSeries(:, N*(i-1)+j), kappa(labels(i, j)), T));
            log_likelihood = log_likelihood + kappa(labels(i,j))*((mean(:, labels(i,j)))')*timeSeries(:, N*(i-1)+j) + log(C_kappa(labels(i, j)));
        end
    end
    
    negative_log_likelihood = -1*log_likelihood;
    outputEnergy = log_prior_spatial + log_prior_label + negative_log_likelihood;
end

