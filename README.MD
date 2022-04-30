# What we have implemented:
We have constructed the whole graph whose minimum cut corresponds to the most optimal alpha-expansion move.
We have also implemented Gibbs sampling and random time-series generation in contiguous and dis-contiguous clusters.
You can change the parameters in randomTester.m and generate a sample input data image with spatially smooth clusters.
However, the edge weights can very well be negative, which restricts us from using the maxflow() function of MatLab.
In fact, no implementation of MatLab code can be found online (even on GitHub) that acomodates -ve edge weights. So what needs to be done in order to see the graph and a random input image is:

# Put the parameters N(size of the NxN image), T(length of the voxel time series), L(number of labels required), value of rho_in(correlation with in the cluster) and rho_ac(correlation across the clusters) in randomTester.m

# Run randomTester.m

# Now call the function AlphaExpansion as follows:
# AlphaExpansion(reshape(timeSeries, [T, N, N]), randi([1, L],N,N), alpha, (1:L), beta_s, beta_l, N, means, kappa, T)
# Here timeSeries is the output of randomTester script. Put the same value of T, N, L in this function call as well. Assign alpha to be a label in the set {1,2,....,L}. beta_s and beta_l are prior parameters(say 0.1 and 0.1 respectively). means is a TxN matrix where each column is unit norm. kappa is a Lx1 column vector of concentration parameters. Design them first according to you and call the function. You will a graph whose mincut corresponds to the min energy of the posterior function.
N = 60
T = 20
L = 5
m = rand([T, L]);
n = sqrt(sum(m.*m));
n=repmat(n, 20, 1);
m = m./n;
kappa = ones(L,1);

example: AlphaExpansion(reshape(timeSeries, [20,60,60]), randi([1, 5],60,60), 1,(1:5), 0.1, 0.1, 60, m, ones(5,1), 20);
# Possible Alternative:
We provide a link to a github project that implemenmts Boykov's original paper.
The MatLab code does not implement the prior on number of image labels, but it does find the minimum cut on a graph with -ve edge weights.
Finally, it also has an example image over which to run the code to test its clutering accuracy.
https://github.com/shaibagon/GCMex

# We wrote the algorithm in fMRI_parcellate.m. Don't run it. The structure of the algorithm is as explained earlier. Had the min cut algorithm been available in MATLAB we would had made a call to it and got the optimal answer. 