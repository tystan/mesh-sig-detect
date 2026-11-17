

# mesh-sig-detect

Time to signal analysis for spontaneous report data on pelvic mesh medical device


![](https://github.com/tystan/mesh-sig-detect/blob/main/fig/sim_multi-grps_multi-test_sig_detect_over_time_thresh-0.05.png)


## Data source

The data is thanks to [curtis-murray](https://github.com/curtis-murray) at his [MedicalDevicesNLP](https://github.com/curtis-murray/MedicalDevicesNLP) repo 


## Repository file structure

* [dat/](https://github.com/tystan/mesh-sig-detect/tree/main/dat) contains the raw and analysis ready data
* [out/](https://github.com/tystan/mesh-sig-detect/tree/main/out) contains processing and analysis output
* [r/](https://github.com/tystan/mesh-sig-detect/tree/main/r) contains `.R` scripts that are prefixed with numbers for rough ordering of steps in processing/analysis/reporting etc
* The `.qmd` files in the top level directory are the [Quarto](https://quarto.org/) (Rmarkdown/Jupiter Notebook-esque documents) that generate the corresponding summary `*.pdf` documents in the top level directory


# Excerpts from Supplementary Material

## Section SM2: Signal detection

Signal detection of disproportionate adverse events (AEs) will often have tabulated count data accumulated over time. The data at time point $t$ can be summarised as below:

|         | AE(s) $\in Y$| AE(s) $\in \bar{Y}$ |
|:--------|---------:|----------:|
|Target exposure|       $a_t$|        $b_t$|
|Comparator exposure|       $c_t$|        $d_t$|

where

* AE(s) $Y$ is the set of AEs (or singular AE) of interest, 
* AE(s) $\bar{Y}$ is the complementary set to the AEs of interest, 
* *Target exposure* is the medical device(s) of interest,  
* *Comparator exposure* is the medical devices to which the *Target exposure* is being compared, and 
* $a_t$, $b_t$, $c_t$ and $d_t$ (all $\in \mathbb{Z}^{+}$) are the respective counts of AEs recorded up until (i.e., cumulative) time $t$.

In the motivating example of the pelvic mesh device, the contingency table can be written more specifically as


|         | Pain AEs | Not pain AEs  |
|:--------|---------:|----------:|
|Pelvic mesh|       $a_t$|        $b_t$|
|Comparator exposure|       $c_t$|        $d_t$|


where 

* *AEs pain* is the count of AEs that contain "pain" themes greater or equal to some pre-specified threshold $p_t \in (0,1)$ as estimated by the hSBM (that is, $P(\texttt{topic == "pain"} | \texttt{Level, Doc}) \geq p_t$), and
* *Comparator exposure* can be any relevant set of medical devices to compare the pelvic mesh to (e.g., hernia mesh or all other mesh devices or all other devices).

### Signal detection over time

We will consider the three signal detection statistics below:

* Proportional reporting ratio (PRR),
* Bayesian Confidence Propagation Neural Network Information Component (BCPNN IC with MCMC CIs), and
* the maxSPRT statistic


As signal detection is being undertaken repeatedly as data are being accumulated, alpha spending needs to be considered. The below table classifies the aforementioned signal detection methods by their null hypothesis as well as how they control for the family-wise error rate (FWER).


| Method | Null hypothesis        |  FWER adjustment |
|:--------|---------:|----------:|
| PRR | Ratio of pain AEs to all AEs in target and comparator groups has a ratio of 1 | PRR with $\alpha$-spending scheme | 
|BCPNN IC |Independence of pain AEs and target group (based on marginal counts) |        IC with $\alpha$-spending scheme| 
| maxSPRT (binary, group sequential)  | Log-likelihood ratio test statistic below null derived threshold | (built in) |


We will demonstrate how the group sequential binary maxSPRT, as described in previous work, is equivalent to a natively FWER-controlled PRR method of signal detection using the log-likelihood ratio test for significance.



Methods used included the Bayesian Confidence Propagation Neural Network (BCPNN) and the maximised Sequential Probability Ratio Test (maxSPRT) which accounted for multiple testing through alpha spending. The BCPNN was used with and without adjusting for multiple testing. The test statistic for BCPNN is the information criterion (IC) which represents the log2 of the ratio of observed to expected adverse events (0 under the null hypothesis of no association between the topic and pelvic mesh).

maxSPRT was developed for near continuous sequential monitoring (called “group sequential” when monitored at discrete time points or after set accumulations of events), maintaining the correct overall alpha level. The test statistic is based on the maximized (log-)likelihood ratio statistic which uses the observed and expected (under the null hypothesis) reporting ratio assuming binomial adverse event accumulation. The critical value of the likelihood ratio statistic is determined by the 100(1-α)% quantile of possible likelihood ratio statistics from binomial adverse event counts under the null hypothesis over the entire group sequential follow-up.






#### Proportional reporting ratio (PRR) 


The PRR estimate is calculated  

```math
\widehat{\text{PRR}}_t = \frac{\frac{a_t}{a_t + b_t}}{\frac{c_t}{c_t+d_t}}.
```

In the context of signal detection, an elevated proportional reporting ratio is of concern. Therefore the one-sided hypothesis test $H_0: \text{PRR} \leq 1$ (proportional reporting of the target is less than the comparator) is used and is not rejected until

$$
\widehat{\text{PRR}}_t {\times} \exp \left\{ - Z_{\alpha}^{*}  \sqrt{ \frac{1}{a_t} + \frac{1}{a_t + b_t} + \frac{1}{c_t} + \frac{1}{c_t+d_t} } \right\} > 1
$$

at the $\alpha$ level where $Z_{\alpha}^{*}$ is the $(1-\alpha)^{\text{th}}$ quantile of the standard normal distribution. The above threshold is equivalent to the lower bound of the approximate $100(1-2\alpha)$\% confidence interval for a standard two-sided hypothesis test.





####  Bayesian Confidence Propagation Neural Network (BCPNN) Information Component (IC)

The Information Component (IC) statistic is an estimate of the observed-to-expected ratio of the number of target exposure AEs of interest on the log${}_2$-scale under independence between the target exposure and AEs of interest based on information theory [(Bate et al., 1998)](https://doi.org/10.1007/s002280050466)
$$
\text{IC}_{XY} = \log_2 \frac{P_{X, Y}(a_t + b_t,a_t + c_t)}{P_X(a_t + b_t) P_Y(a_t + c_t)}
$$
where $P_X(X=x)$ denotes the marginal probability of an observed count $x$ for the target exposure, $P_Y(Y=y)$ denotes the marginal probability of an observed count $y$ for the AE of interest, and $P_{X,Y}(X=x, Y=y)$ denotes the joint probability.

The BCPNN IC of [Noren et al. (2006)](https://doi.org/10.1002/sim.2473) uses a Bayesian inference based *maximum a posteriori* (m.a.p.) central estimate of the IC,  
$$
\widehat{\text{IC}}_t = \log_2 \frac{\text{E}\left[ \hat{p}_a \right]}{\text{E}\left[ \hat{p}_a + \hat{p}_b\right]\text{E}\left[ \hat{p}_a + \hat{p}_c\right]}
$$

where $p_a$, $p_b$ and $p_c$ are the (assumed constant over time) underlying probabilities of the multinomial-distributed observed events $a_t$, $b_t$ and $c_t$, respectively ($p_d$ corresponding to the count $d_t$ also included). The underlying probabilities are modelled using Dirichlet priors resulting in a Dirichlet posterior distribution. The one-sided null hypothesis of the joint probability target exposure and AEs of interest is equal or less than the marginal products ($H_0: \text{IC}_t \leq 0$) can be rejected when the $\alpha$ quantile of the Markov Chain Monte Carlo (MCMC) empirical distribution is greater than 0. Similarly to the rejection rule for the PRR, this threshold corresponds to the lower bound of the $100(1-2\alpha)$\% equal-tailed credible region in a two-sided hypothesis test.



#### maxSPRT

[Kulldorff et al. (2011)](https://doi.org/10.1080/07474946.2011.539924) outlined that the relative risk (RR) at a given point-in-time for accumulated binary data (that is, "success"/"failure" events or AE of interest or not) of a target group relative to a comparator has the maximum likelihood estimate of 

$$
\widehat{\text{RR}} = z \frac{C_n}{n - C_n}
$$

where 

* $z$ is the ratio of the total AEs for the comparator to the total AEs for the target,
* $C_n$ is the count of target exposure AEs in $X$,
* $n$ is the count of all AEs in $X$ (target and comparator exposure), and
* $n - C_n$ is therefore the count of comparator exposure AEs in $X$.

In the context of our data, the values $z$, $C_n$ and $n$ are the quantities $\frac{c_t + d_t}{a_t + b_t}$, $a_t$ and $a_t + c_t$, respectively, at time $t$.

Therefore the RR maximum likelihood estimate at time $t$ can be re-written
$$
\begin{aligned}
\widehat{\text{RR}}_t 
&= 
\frac{c_t + d_t}{a_t + b_t} \times \frac{a_t}{c_t}
\\
&= 
\frac{\frac{1}{a_t + b_t}}{\frac{1}{c_t + d_t}} \times \frac{a_t}{c_t}
\\
&= 
\frac{\frac{a_t}{a_t + b_t}}{\frac{c_t}{c_t + d_t}} 
\\
\end{aligned}
$$
which is the PRR estimate at time $t$ as before.

The (maximised) log-likelihood ratio statistic of $\widehat{\text{PRR}}_t$ (equivalently, $\widehat{\text{RR}}_t$) can be determined calculated as

$$
\text{LLR}_t = 
a_t  \ln \left(\frac{a_t}{a_t + c_t}\right) + 
c_t  \ln \left(\frac{c_t}{a_t + c_t} \right) - 
a_t  \ln\left(\frac{a_t + b_t}{a_t + b_t +c_t + d_t}\right) - 
c_t \ln\left(\frac{c_t + d_t}{a_t + b_t + c_t + d_t}\right)
$$

The maxSPRT test is considered significant when $LLR_t$ is greater than the pre-computed critical value which is the $100(1-\alpha)$\% percentile of the  $LLR_t$ values generated under the null hypothesis $\text{RR}_t = 1$ for group sequential looks at the data $t = t_1, t_2, \ldots, t_k$. The CV can either be computed using the 95th percentile of $LLR_t$ values with the exact joint binomial probabilities over $k$ looks of the data accumulation under the null hypothesis, or by MCMC sampling of binomial event accumulation (and associated $LLR_t$ values) to approximate the $LLR$ distribution when the exact CV computation is computationally intractable.\footnote{The \texttt{Sequential} R package exact CV function suggests "not using values greater than 1000" in regards to the total events to potentially be observed over the follow-up.}






