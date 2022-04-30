function [G] = AlphaExpansion(timeSeries, initialLabels, alpha, L_in, beta_s, beta_l, N, means, kappa, T)
    % ALPHAEXPANSION - uses Boykov's method to minimize energy
    % Uses a slightly different graph to also incorporate label energy
    numOfLabels = size(L_in);
    numOfLabels = numOfLabels(2);
    edges_1 = [];
    edges_2 = [];
    edge_weights = [];
    sum_edges = zeros(1, 100);
    C_kappa = kappa.^(T/2 - 1);
    C_kappa = C_kappa./(besseli(T/2 - 1,kappa));
    C_kappa = C_kappa/((2*pi)^(T/2));
    for i=1:N
        for j=1:N
            edges_2 = [edges_2, "p_"+string(i)+"_"+string(j)];
            if (initialLabels(i, j) == alpha)
                edges_1 = [edges_1, "alpha"];
            else
                edges_1 = [edges_1, "l_"+string(initialLabels(i, j))];
            end
            edge_weight = - (log(C_kappa(alpha)) + (kappa(alpha)*((means(:, alpha))')*timeSeries(:, i, j)));
            %prod = (means(:, alpha)')*timeSeries(:, i, j);
            %edge_weight = -log(VMFMeanDirDensity(prod, kappa(alpha), T));
            edge_weights = [edge_weights, edge_weight];
            sum_edges(initialLabels(i, j)) = sum_edges(initialLabels(i, j)) + edge_weight;
            edges_2 = [edges_2, "p_"+string(i)+"_"+string(j)];
            edges_1 = [edges_1, "alpha_bar"];
            edge_weight = - (log(C_kappa(initialLabels(i, j))) + (kappa(initialLabels(i, j))*((means(:, initialLabels(i, j)))')*timeSeries(:, i, j)));
            %prod = (means(:, alpha)')*timeSeries(:, i, j);
            %edge_weight = -log(VMFMeanDirDensity(prod, kappa(initialLabels(i, j)), T));
            edge_weights = [edge_weights, edge_weight];
        end
    end
    for i=1:numOfLabels
        thisLabel = L_in(i);
        if (thisLabel == alpha)
            continue;
        end
        edges_1 = [edges_1, "alpha"];
        edges_2 = [edges_2, "l_"+string(thisLabel)];
        edge_weights = [edge_weights, sum_edges(thisLabel) - beta_l];
    end
    for i=1:N
        for j=1:N
            if (initialLabels(i, j) ~= initialLabels(mod(i, N)+1, j))
                edges_1 = [edges_1, "p_"+string(i)+"_"+string(j)];
                edges_2 = [edges_2, "a_"+string(i)+"_"+string(j)+"_d"];
                edge_weights = [edge_weights, beta_s*(alpha~=initialLabels(i, j))];
                edges_1 = [edges_1, "p_"+string(mod(i, N)+1)+"_"+string(j)];
                edges_2 = [edges_2, "a_"+string(i)+"_"+string(j)+"_d"];
                edge_weights = [edge_weights, beta_s*(alpha~=initialLabels(mod(i, N)+1, j))];
                edges_1 = [edges_1, "alpha_bar"];
                edges_2 = [edges_2, "a_"+string(i)+"_"+string(j)+"_d"];
                edge_weights = [edge_weights, beta_s*(initialLabels(mod(i, N)+1, j)~=initialLabels(i, j))];
            else
                edges_1 = [edges_1, "p_"+string(i)+"_"+string(j)];
                edges_2 = [edges_2, "p_"+string(mod(i, N)+1)+"_"+string(j)];
                edge_weights = [edge_weights, beta_s*(alpha~=initialLabels(i, j))];
            end
        end
    end
    for i=1:N
        for j=1:N
            if (initialLabels(i, j) ~= initialLabels(i, mod(j, N)+1))
                edges_1 = [edges_1, "p_"+string(i)+"_"+string(j)];
                edges_2 = [edges_2, "a_"+string(i)+"_"+string(j)+"_r"];
                edge_weights = [edge_weights, beta_s*(alpha~=initialLabels(i, j))];
                edges_1 = [edges_1, "p_"+string(i)+"_"+string(mod(j, N)+1)];
                edges_2 = [edges_2, "a_"+string(i)+"_"+string(j)+"_r"];
                edge_weights = [edge_weights, beta_s*(alpha~=initialLabels(i, mod(j, N)+1))];
                edges_1 = [edges_1, "alpha_bar"];
                edges_2 = [edges_2, "a_"+string(i)+"_"+string(j)+"_r"];
                edge_weights = [edge_weights, beta_s*(initialLabels(i, mod(j, N)+1)~=initialLabels(i, j))];
            else
                edges_1 = [edges_1, "p_"+string(i)+"_"+string(j)];
                edges_2 = [edges_2, "p_"+string(i)+"_"+string(mod(j, N)+1)];
                edge_weights = [edge_weights, beta_s*(alpha~=initialLabels(i, j))];
            end
        end
    end
    %figure();
    G = graph(edges_1, edges_2, edge_weights);
    plot(G);
%     plot(G, 'EdgeLabel', G.Edges.Weight);
%     figure();
%     [~, ~, cs, ct] = maxflow(G, "alpha", "alpha_bar");
%     for i=1:numOfLabels
%         thisLabel = L_in(i);
%         if (thisLabel == alpha)
%             L_out = [L_out, alpha];
%             continue;
%         end
%         if (any(ct(:)=="l_"+string(thisLabel)))
%             indexMatrix = (thisLabel==initialLabels);
%             finalLabels(indexMatrix) = alpha;
%         else
%             L_out = [L_out, thisLabel];
%         end
%     end
%     for i=1:N
%         for j=1:N
%             if (any(cs(:)=="p_"+string(i)+"_"+string(j)))
%                 finalLabels(i, j) = initialLabels(i, j);
%             else
%                 finalLabels(i, j) = alpha;
%             end
%         end
%     end
end