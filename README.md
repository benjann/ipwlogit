# ipwlogit
Stata module to fit marginal logistic regression by inverse probability weighting

`ipwlogit` fits marginal logistic regression of a binary dependent variable on
a treatment variable, possibly adjusting for control variables by inverse
probability weighting (IPW). The resulting estimate can be interpreted as a
marginal (log) odds ratio of a positive outcome. The treatment variable can be
categorical, continuous, or discrete.

---

Installation from GitHub:

    . net install ipwlogit, replace from(https://raw.githubusercontent.com/benjann/ipwlogit/main/)

`ipwlogit` with option `psmethod(gologit)` requires
[`gologit2`](https://www3.nd.edu/~rwilliam/gologit2/) by Richard Williams. Type

    . ssc install gologit2

to install `gologit2`.

---

Main changes:

    15aug2022 (version 1.0.0):
    - ipwlogit released on GitHub