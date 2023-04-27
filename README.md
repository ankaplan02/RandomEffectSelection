# Random Effect Selection
Chen and Dunson (2003) provided a means of random effect variable selection via a version of Cholesky Decomposition of the random effects in a linear mixed modeling framework. The variable selection priors are placed on the diagonal entries, assumed to be zero-inflated truncated normal distribution on the non-negative real line. 

In certain research contexts it may be critical in assessing whether the effects of certain sample characteristics on the outcome vary between groups. Moreover, it may be important to assess how important these between-group differences really are. 
Chen and Dunson (2003) illustrate a novel extension of the linear mixed effects model for continuous outcomes that evaluates the probability of a random effect variance is zero. 

Let $y_i$ denote the outcome of group $i$ such that this vector $y_i = (y_{i1},...,y_{in_{i}})^T$, with covariate matrix $X_i = (x_{i1}^T,..., x_{in_{i}}^T)^T$ where each $x_{ij}$ is a row vector of length $q_1$. Follow this with a group-specific vector of predictors $Z_{i} = (z_{i1}^{T},...,z_{in_{i}}^T)^T$, this matrix is $n_i \times q_2$ where $q_2$ is the number of random effects. Then the following linear mixed regression can be used: 

$$y_i = X_i \alpha + Z_i \beta_i + \epsilon_i$$

Typically, $\alpha$ is the $q_1 \times 1$ vector of unknown study-population parameters while $\beta_i \sim N(0,D)$, and the residual vector $\epsilon_{i} \sim N(0, \sigma^2 \times I_{n_i})$ where $I_{r}$ refers to the identity matrix of size $r \times r$. $D$ is $q_2 \times q_2$. 

The reparameterization specifies that 

$$y_i = X_i \alpha + Z_i \Lambda \Gamma b_i + \epsilon_i$$

where $\Lambda \times \Gamma \times b_i = \beta_i$, and $\Lambda$ is a diagonal matrix with entries $\lambda = (\lambda_1, \lambda_2,...,\lambda_{q_2})^T$ and $\Gamma$ is lower triangular with 
diagonal entries equal to 1, and $q_{2}\times (q_{2}-1)/2$ lower triangular entries $\gamma_{ml}$ where $m = 2,...,q_{2}$ and $l = 1,...,q_{2}-1$. This reformulation has the key equality that 
$D = \Lambda \Gamma \Gamma^T \Lambda$, so the on-outcome scale variability measured by the diagonal entries in $D$ and the correlations between the random effects by the off-diagonal entries in $D$ can be obtained. 

The prior specifications are primarily conjugate priors and are as follows:

$$\alpha \sim N(\alpha_0, A_0)$$
$$b_i \sim N(0, I_{q_2})$$
$$\sigma^{-2} \sim Gamma(c_0, d_0)$$
and $c_0$ and $d_0$ are specified to be 0.05 by default but you may change it in the function call; 
$$\gamma \sim N(\gamma_0, R_0)$$
where these are the lower diagonal entries in $\Gamma$ - in the provided code the elements of $\gamma$ are in the following order $(\gamma_{21}, \gamma_{31}, \gamma_{32}, \gamma_{41}, \gamma_{42},\gamma_{43},...)$ and so on. This is to retain the structure of the lower triangular matrix how Chen and Dunson specified it to be; 
$$\lambda_l \sim ZI-N^{+}(p_{l0}, m_{l0}, s_{l0}^2)$$
where the $p_{l0}$ represents the *prior* probability that $\lambda_l$ is 0, $m_{l0}$ and $s_{l0}$ are the prior mean and standard deviation for $\lambda_l$. In the article, the authors provide various specifications for these parameters. The notation $ZI-N^{+}$ denotes the Zero-Inflated Truncated Normal distribution restricted to be on the non-negative real line.  

The supplied code is all in R but future subroutines done in C++ may be of huge benefit by speeding up vectorization and for-loops to analyze larger data sets. 
