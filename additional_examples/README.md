# Examples of computing divergent master integrals

In the directory `additional_examples`, we provide proof-of-concept examples to show how to use semi-definite programming to compute master integrals that are either ultraviolet (UV) or infrared (IR) divergent. 


## Computing UV divergent master integrals

If the concerned master integrals are known to have potential UV divergences only (it is the case, for example, when all propagators have non-zero masses, and in the Euclidean region all other divergences are forbidden by kinematics), then their computation can be handled without introducing too much further trick. The basic idea is to first lower the space-time dimension $d$ by an even integer to make the UV divergences disappear, and then use **dimensional-shifting relations** to move the space-time dimension upwards.

As an example, we consider the three-loop unequal-mass banana integral family already treated in [Feynman Integrals from Positivity Constraints](https://arxiv.org/pdf/2303.15624), but at a space-time dimension of $d=4$ where some master integrals have UV divergences.

Compile and run `additional_examples/uvexample.cpp` by the command `make uv`, and the program will first output a dimensional-shifting matrix, which transforms $(d-2)$-dimensional master integrals into $d$-dimensional master integrals through a left multiplication. Then the program solves the master integrals at space-time dimension $d=0$. One needs to manually left multiply the dimensional-shifting matrix on the result twice to produce the desired result.

## Computing IR divergent master integrals

Master integrals with potential IR divergences are harder to deal with. Simply raising the space-time dimension does not always work, as it may introduce potential UV divergences. In this case, we find that IBP relations alone are actually not enough to fully determine the values of master integrals.

The strategy to deal with IR divergences is to combine the semi-definite programming technique with differential equations. The computation is decomposed into the following steps:

* First, introduce an auxiliary mass $x>0$ for all (or some) massless propagators, such that all the master integrals are IR finite. Notice that this process may increase the number of master integrals.

    It is not hard to see that introducing $x$ does not break Euclidean condition, and that IR divergences in the original master integrals are manifest through a singularity at $x=0$ in the complex plane of $x$. 

* Then we build the differential equation of master integrals with respect to kinematic variable $x$. It is well known that $x=0$ is at most a regular singularity of the equation.

    There are well-developed algorithms to compute the general solutions of a system of linear differential equations around its regular singularity by power series method. We give an implementation [here](https://github.com/zml72062/linear-ode), which has the advantage of being able to solve differential systems depending on the dimensional regularizer $\epsilon$. We use it to generate the general solutions of the above differential equation for master integrals. If there are $N$ master integrals after introducing $x$, then the general solutions constitute a $N\times N$ non-singular matrix $S(x,\epsilon)$.

* Now it suffices to determine the arbitrary constants $c(\epsilon)$ in the general solutions. We achieve this by matching the result of semi-definite programming and the solution of differential equation. Assume that $x=x_0>0$ lies within the radius of convergence of the series solutions. We compute the $\epsilon$-expansion of $c(\epsilon)$ in two steps,

    * we first plug $x=x_0$ into the general solution $S(x,\epsilon)$ of differential equation, and compute its matrix inverse $S^{-1}(x_0,\epsilon)$ in terms of a Laurent expansion in $\epsilon$ 
    
    * we then use semi-definite programming to compute the $\epsilon$-expansion of master integrals at $x=x_0$, which should now succeed since the master integrals are convergent for any $x=x_0>0$

    Finally, we left multiply $S^{-1}(x_0,\epsilon)$ on the column vector consisting of master integral values at $x=x_0$. We expect the result to be $c(\epsilon)$ and independent of the choice of $x_0$.

* Eventually, we take the limit $x\rightarrow 0$ in $S(x,\epsilon)\cdot c(\epsilon)$, and set to zero all terms proportional to $x^{a+b\epsilon}$ with $b\ne 0$. This would be the desired result. 

As an example, we still consider the three-loop unequal-mass banana integral family at $d=2$, but with one mass $m_1$ set to zero. In this case, all master integrals are free from UV divergences, but some have IR divergences as loop momentum $\ell_1$ is small. The first step to solve the master integrals is to give the propagator $\ell_1^2$ a non-zero mass $x$, namely $\ell_1^2\rightarrow \ell_1^2-x$.

First compile `additional_examples/irexample.cpp` by the command `make ir`. Running the program with no command-line argument will print the coefficient matrix of the differential equation with respect to $x$. One needs to provide this matrix to the linear differential equation [solver](https://github.com/zml72062/linear-ode) to get the general solution $S(x,\epsilon)$. This correponds to the second step of the solving procedure.

One then runs the program again, but with a command-line argument selected from `0.005`, `0.010`, `0.015` or `0.020`. This will use semi-definite programming to evaluate the $\epsilon$-expansion of master integrals at $x=0.005$, $x=0.010$, $x=0.015$ or $x=0.020$. 

After the above results are obtained, one can use a symbolic manipulator to find out $S^{-1}(x,\epsilon)$ at $x=0.005$, $x=0.010$, $x=0.015$ or $x=0.020$, and then determine $c(\epsilon)$ and the final result. The procedure of plugging in numeric values of $x$, computing inverse and getting the final result needs to be carried out manually, and we save a copy of intermediate results in `additional_examples/cache/ir`. The definitions of files are described below.

|File | Description |
|:---:|:---|
|`diffeq`| Coefficient matrix of differential equation. |
|`radius_convergence`| Radius of convergence of power series solution around $x=0$.|
|`diffeq_sols`| General solution $S(x,\epsilon)$ at four non-zero $x$ points.|
|`inv_diffeq_sols`| Inverse of general solution $S^{-1}(x,\epsilon)$ at four non-zero $x$ points.|
|`diffeq_sol0`| General solution $S(x,\epsilon)$ at $x\rightarrow 0$. |
|`sdp_results`| Semi-definite programming results at four non-zero $x$ points.|
|`raw_result`| The final result, which contains spurious poles in $\epsilon$ (such as `1e-80/eps^4`) due to numerical error. |
|`result`| The final result after removing spurious poles and truncating to desired order. |


