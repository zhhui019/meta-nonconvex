Meta-Analysis Based on Nonconvex Regularization

*****************************************************************************************************************
* author: Zhang Hui  											*
* email:    zhanghui.nwu@foxmail.com
* Last revision: 01-March-2019                                  								*
*****************************************************************************************************************
* RECOMMENDATIONS:                                   	*
* This toolbox is designed to work with Matlab 2018a  *
*********************************************************

------------------------------------------------------------------------------------------------------------------------------------------------
DESCRIPTION:
This toolbox provides an efficient way to Solve the meta-analysis method 
based on nonconvex regularization problem with a single tuning parameter.
Jointly fit a logistic regression model with a nonconvex penalty over multiple datasets.
It enables heterogeneous variable selections in different datasets. 
The function minimizes \eqn{-logLik + lambda * p(beta)}, where \eqn{-logLik} is the negative of the total
log-Likelihood from all datasets, \eqn{lambda} is a single tuning parameter and \eqn{p(beta)} is a specific nonconvex penalty
function enabling heterogeneous selections of variables in different datasets. For more details of the penalty
function, see the reference below.
------------------------------------------------------------------------------------------------------------------------------------------------
SPECIFICATIONS for using meta_nonconvex

One demo file 'demo_meta_nonconvex.m' is provided to show how to use the 'meta_nonconvex' algorithm.

The main function is 'meta_nonconvex.m'.

'Logistic_nonconvex_func.m' is to solve the subproblem (logistic regression model with nonconvex penalty).

------------------------------------------------------------------------------------------------------------------------------------------------
RELATED REFERENCES:

[1] Hui Zhang, Shou-Jiang Li, Hai Zhang, Zi-Yi Yang, Yan-Qiong Ren, Liang-Yong Xia and Yong Liang,
    Meta-Analysis Based on Nonconvex Regularization.
[2] Hui, Zhang, Hai, et al. Approximate Message Passing Algorithm for Nonconvex Regularization[J]. IEEE Access, 2019.
[3] 张会, 张海. 基于AMP的L_(1/2)正则化方法[J]. 中国科学:信息科学, 2017(01):62-76.

------------------------------------------------------------------------------------------------------------------------------------------------
