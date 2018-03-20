Small introduction: 

**contrast_experiments** is a collection of experiments devoted to contrast estimation in FitzHugh-Nagumo model (that's how we get the results from the preprint).

**dimensionality_test_experiments** are more recent, and there are few pretty cool functions, in particular:
+ *construct_test* (implementation of statistical tests from paper by Jean Jacod and Mark Podolskij "A test for the rank of the volatility process: the random perturbation approach"). Returns the value of the dimensionality estimator $\hat{R}$ and the value of the test statistics $V$. 
+ *hawkes_approximation* (implementation of the stochastic diffusion which approximates the Hawkes process, from paper by Eva Löcherbach and Susanne Ditlevsen "Multi-class oscillating systems of interacting neurons)
+ *FHN.simulate* (simulations of FitzHugh-Nagumo model with higher-order scheme)

**example.cpp** few matrix algebra functions, implemented in Rcpp/RcppEigen, which are used to speed up the code (which is pretty costly in terms of computation!)

If you have noticed any inaccuracies - you can contact me (up-to-date contact information is available on my [website](http://amelnykova.com)). 
