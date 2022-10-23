# ipwlogit
Stata module to fit marginal logistic regression by inverse probability weighting

`ipwlogit` fits marginal logistic regression of a binary dependent variable on
a treatment variable, possibly adjusting for control variables by inverse
probability weighting (IPW). The resulting estimate can be interpreted as a
marginal (log) odds ratio of a positive outcome. The treatment variable can be
categorical, continuous, or discrete.

---

Installation from GitHub:

    . net from https://raw.githubusercontent.com/benjann/ipwlogit/main/
    . net install ipwlogit, replace

`ipwlogit` with option `psmethod(gologit)` requires
[`gologit2`](https://www3.nd.edu/~rwilliam/gologit2/) by Richard Williams. Type

    . ssc install gologit2

to install `gologit2`.

---

Main changes:

    20oct2022 (version 1.0.7)
    - truncate() is now applied to the overall distribution of stabilized weights
      across all treatment levels, i.e. truncate() is no longer applied to each
      treatment levels individually; non-stabilized weights are computed after
      applying truncate()
    - e(ipwstats) renamed to e(ipw)
    - marginal treatment probabilities now returned in e(prop)

    20oct2022 (version 1.0.6)
    - option truncate() added

    04sep2022 (version 1.0.5)
    - tvar may now contain polynomials; parsing of varlist improved
    - for categorical tvar, IPWs will now always be based on the observed levels,
      not the levels specified in tvar
    - options -rifgenerate()- and -ifscaling()- added
    - e(sum_w) added to returns

    01sep2022 (version 1.0.4)
    - now using Sturges' rule to determine the number of bins used to categorize
      a continuous treatment

    17aug2022 (version 1.0.3)
    - now requires Stata 14 or newer
    - ttype "categorical" renamed to "factor"
    - some adjustments to header display

    16aug2022 (version 1.0.2)
    - option nodots added

    15aug2022 (version 1.0.1)
    - fweights now supported
    
    15aug2022 (version 1.0.0):
    - ipwlogit released on GitHub