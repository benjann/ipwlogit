*! version 1.0.4  01sep2022  Ben Jann

program ipwlogit, eclass properties(svyb svyj mi)
    version 14
    if replay() {
        Display `0'
        exit
    }
    local version : di "version " string(_caller()) ":"
    _parse_opts `0' // returns 00
    `version' _vce_parserun ipwlogit, mark(CLuster) : `00'
    if "`s(exit)'" != "" {
        ereturn local cmdline `"ipwlogit `0'"'
        exit
    }
    _ipwlogit `00' // returns diopts
    eret local cmdline `"ipwlogit `0'"'
    Display, `diopts'
    if `"`e(generate)'`e(tgenerate)'`e(ifgenerate)'"'!="" {
        describe `e(generate)' `e(tgenerate)' `e(ifgenerate)'
    }
end

program _parse_opts
    _parse comma lhs 0 : 0
    syntax [, or eform(passthru) vce(passthru) NOVCEADJust ///
        GENerate(passthru) TGENerate(passthru) IFGENerate(passthru) * ]
    local options `vce' `options'
    // replace option -or- by eform()
    if "`or'"!="" {
        if `"`eform'"'=="" local eform eform(Odds ratio)
        local options `eform' `options'
    }
    // bootstrap/jackknife
    _parse_opts_vce, `vce' // returns repl
    if `repl' {
        if `"`generate'"'!="" {
            di as err "{bf:generate()} not allowed with replication based VCE"
            exit 198
        }
        if `"`tgenerate'"'!="" {
            di as err "{bf:tgenerate()} not allowed with replication based VCE"
            exit 198
        }
        if `"`ifgenerate'"'!="" {
            di as err "{bf:ifgenerate()} not allowed with replication based VCE"
            exit 198
        }
        local novceadjust novceadjust
    }
    // returns
    c_local 00 `lhs', `generate' `tgenerate' `ifgenerate' `novceadjust' `options'
end

program _parse_opts_vce
    syntax [, vce(str) ]
    gettoken v : vce, parse(", ")
    local repl 0
    if      `"`v'"'==substr("bootstrap",1,max(4,strlen(`"`v'"'))) local repl 1
    else if `"`v'"'==substr("jackknife",1,max(4,strlen(`"`v'"'))) local repl 1
    c_local repl `repl'
end

program Display
    syntax [, or noHEADer noTABle noIPWstats * ]
    if "`or'"!="" local eform eform(Odds Ratio)
    if "`header'"=="" {
        local w1 17
        local c1 49
        local c2 = `c1' + `w1' + 1
        local w2 10
        local headopts head2left(`w1') head2right(`w2')
        if      c(stata_version)<17            local headopts
        else if d(`c(born_date)')<d(13jul2021) local headopts
        _coef_table_header, `headopts'
        di as txt _col(`c1') "Treatment type" _col(`c2') "= " as res %`w2's e(ttype)
        di as txt _col(`c1') "Number of levels" _col(`c2') "= " as res %`w2'.0g e(tk)
        di as txt _col(`c1') "PS method" _col(`c2') "= " as res %`w2's e(psmethod)
        di ""
    }
    if "`table'"=="" {
        eret display, first `eform' `options'
        local xvars `"`e(indepvars)'"'
        if `"`xvars'"'!="" {
            if (strlen(`"`xvars'"')>63) {
                local xvars = substr(`"`xvars'"',1,60) + "..."
            }
            di as txt "(adjusted for " as res `"`xvars'"' as txt ")"
        }
    }
    if "`ipwstats'"=="" {
        di _n as txt "Distribution of IPWs" _c
        matlist e(ipwstats), lines(none) showcoleq(combined) row(level)
    }
end

program _ipwlogit, eclass
    syntax varlist(min=1 numeric fv) [if] [in] [pw fw iw] [, ///
        PSMethod(str) PSOpts(str asis) ///
        NOBINARY /// undocumented; use general code even if treatment is binary
        BINs(numlist int max=1 >=2) DISCRete ASBALanced ///
        GENerate(name) TGENerate(name) IFGENerate(str) replace ///
        vce(passthru) Robust CLuster(passthru) NOVCEADJust ///
        NOIsily eform(passthru) noHEADer noIPWstats noTABle ///
        noCONStant noDOTs * ]
    
    // generate option
    if "`generate'`tgenerate'"!="" & "`replace'"=="" {
        confirm new variable `generate' `tgenerate'
    }
    
    // psmethod option
    _parse_psmethod, `psmethod' psopts(`psopts') // returns psmethod
    
    // collect diopts
    _get_diopts diopts options, `options'
    c_local diopts `eform' `header' `ipwstats' `table' `diopts'

    // maximize options (for outcome model)
    mlopts mlopts, `options'
    
    // mark sample and parse vce option
    marksample touse
    _vce_parse `touse', optlist(Robust) argoptlist(CLuster) old: ///
        [`weight'`exp'], `vce' `robust' `cluster'
    local vce      `"`r(vceopt)'"'
    local clustvar `"`r(cluster)'"'
    
    // prepare weights
    if "`weight'"!="" {
        local wvar = substr(`"`exp'"', 3, .) // strip "= "
        capt confirm var `wvar'
        if _rc==1 exit _rc
        if _rc {
            tempname wvar
            qui gen double `wvar' `exp' if `touse'
        }
        // weights for _pctile, summarize, tabstat
        if      "`weight'"=="pweight" local awgt `"[aw=`wvar']"'
        else if "`weight'"=="iweight" local awgt `"[aw=`wvar']"'
        else                          local awgt `"[`weight'=`wvar']"'
        // weights for PS estimation
        if      "`weight'"=="pweight" local iwgt `"[iw=`wvar']"'
        else                          local iwgt `"[`weight'=`wvar']"'
        // weights for total
        local twgt `"[`weight'=`wvar']"'
    }
    
    // iweights imply novceadjust
    if "`weight'"=="iweight" local novceadjust novceadjust
    local vceadj = "`novceadjust'"==""
    
    // parse variables
    gettoken depvar indepvars : varlist
    if `:list depvar in indepvars' {
        di as err "{it:depvar} must be unique"
        exit 198
    }
    gettoken tvar indepvars : indepvars, bind
    if `:list tvar in indepvars' {
        di as err "{it:tvar} must be unique"
        exit 198
    }
    local indepvars `indepvars'
    _fv_check_depvar `depvar'
    fvexpand `tvar' if `touse'
    local texpand `"`r(varlist)'"'
    if `"`r(fvops)'"'!="" local ttype "factor"
    else                  local ttype "continuous"
    if "`ttype'"=="factor" {
        if "`bins'"!="" {
            di as err "{bf:bins()} not allowed with categorical treatment"
            exit 198
        }
        if "`discrete'"!="" {
            di as err "{bf:discrete} not allowed with categorical treatment"
            exit 198
        }
        if "`tgenerate'"!="" {
            di as err "{bf:tgenerate()} not allowed with categorical treatment"
            exit 198
        }
    }
    else {
        if "`discrete'"!="" {
            if "`bins'"!="" {
                di as err "{bf:bins()} and {bf:discrete} not both allowed"
                exit 198
            }
            local ttype "discrete"
        }
    }
    fvexpand `indepvars' if `touse'
    local xvars `"`r(varlist)'"'
    
    // parse ifgenerate()
    if `"`ifgenerate'"'!="" {
        if strpos(`"`ifgenerate'"',"*") {
            gettoken ifstub rest : ifgenerate, parse("* ")
            if `"`rest'"'!="*" {
                di as err "ifgenerate() invalid; " /*
                    */ "must specify {it:stub}{bf:*} or {it:namelist}"
                exit 198
            }
            confirm name `ifstub'
            local k: list sizeof texpand
            if "`constant'"=="" local k = `k' + 1
            local ifgenerate
            forv i = 1/`k' {
                local ifgenerate `ifgenerate' `ifstub'`i'
            }
        }
        confirm names `ifgenerate'
        if "`replace'"=="" {
            confirm new variable `ifgenerate'
        }
    }
    
    // process outcome variable
    tempvar Y
    qui gen byte `Y' = (`depvar'!=0) if `touse'
    sum `Y' if `touse' `awgt', meanonly
    local N = r(N)
    if `N'==0 {
        error 2000
    }
    if r(min)==r(max) {
        di as err "outcome does not vary"
        di as err "remember:                           0 = negative outcome"
        di as err "          all other nonmissing values = positive outcome"
        exit 2000
    }
    
    // default number of bins
    if !inlist("`ttype'","factor","discrete") {
        if "`bins'"=="" local bins = max(1, ceil(ln(`N')/ln(2)) + 1)
    }
    
    // process treatment variable
    // - categorical
    if "`ttype'"=="factor" {
        _collect_fvinfo `texpand' // returns tname, tk, tlevels, tbase
        local T `tname'
        if `tk'<2 {
            di as err "treatment does not vary"
            exit 2000
        }
        if "`psmethod'"=="" {
            if `tk'==2 local psmethod "logit"
            else       local psmethod "mlogit"
        }
    }
    // - continuous
    else {
        local tname "`texpand'"
        if "`psmethod'"=="" local psmethod "cologit"
        // categorize
        tempvar T
        tempname AT
        if "`ttype'"=="discrete" {
            mata: _mkgroups("`tname'", "`T'", "`touse'", "tlevels", "`AT'")
        }
        else {
            qui gen byte `T' = 1 if `touse'
            _pctile `tname' if `touse' `awgt', nquantiles(`bins')
            local tk = `bins' - 1
            local tlevels 1
            local l 1
            forv i = 1/`tk' {
                if r(r`i')==r(r`=`i'+1') continue // skip ties
                if `i'==`tk' { // last cut
                    capt assert (`tname'<=r(r`i')) if `touse'
                    if _rc==1 exit _rc
                    if _rc==0 {
                        mat `AT' = nullmat(`AT'), r(r`i')
                        continue, break // no obs above last cut
                    }
                }
                local ++l
                qui replace `T' = `l' if `tname'>r(r`i') & `touse'
                mat `AT' = nullmat(`AT'), r(r`i')
                local tlevels `tlevels' `l'
                if `i'==`tk' {
                    mat `AT' = nullmat(`AT'), . // upper bound if last bin = .
                }
            }
            mat coln `AT' = `tlevels'
        }
        local tk : list sizeof tlevels
        if `tk'<2 {
            di as err "(coarsened) treatment does not vary"
            exit 2000
        }
        local tbase 1
    }
    
    // compute IPWs
    if "`noisily'"=="" {
        local qui quietly
        if "`dots'"=="" di as txt "(estimating balancing weights " _c
    }
    tempvar ipw
    qui gen double `ipw' = .
    // - logit
    if "`psmethod'"=="logit" {
        tempvar tmpT q
        if `tk'==2 & "`nobinary'"=="" {     // case 1: single logit
            if "`noisily'`dots'"=="" di as txt "..." _c
            qui gen byte `tmpT' = `tname'!=`tbase' if `touse'
            `qui' logit `tmpT' `xvars' `iwgt' if `touse', `psopts'
            qui predict double `q' if `touse', pr
            qui replace `q' = 1 - `q' if !`tmpT' & `touse'
            if "`asbalanced'"!="" {
                qui replace `ipw' = 1 / `q' if `touse'
            }
            else {
                su `tmpT' `awgt' if `touse', meanonly
                qui replace `ipw' = cond(`tmpT', r(mean), 1-r(mean)) / `q' if `touse'
            }
            if `vceadj' {
                tempname tV sc
                _get_V `tV'
                qui predict double `sc' if e(sample), score
            }
        }
        else {                      // case 2: individual logit for each level
            qui gen byte `tmpT' = .
            qui gen byte `q' = .
            local tV
            local sc
            tempvar ql
            foreach l of local tlevels {
                if "`noisily'`dots'"=="" di "." _c
                else               di _n as res "==> level `l'"
                qui replace `tmpT' = `T'==`l' & `touse'
                `qui' logit `tmpT' `xvars' `iwgt' if `touse', `psopts'
                qui predict double `ql' if `tmpT', pr
                qui replace `q' = `ql' if `tmpT'
                drop `ql'
                if "`asbalanced'"!="" {
                    qui replace `ipw' = 1 / `q' if `tmpT'
                }
                else {
                    su `tmpT' `awgt' if `touse', meanonly
                    qui replace `ipw' = r(mean) / `q' if `tmpT'
                }
                if `vceadj' {
                    tempname tVl scl
                    local tV `tV' `tVl'
                    local sc `sc' `scl'
                    _get_V `tVl'
                    qui predict double `scl' if e(sample), score
                }
            }
        }
    }
    // - mlogit
    else if "`psmethod'"=="mlogit" {
        if "`noisily'`dots'"=="" di "..." _c
        `qui' mlogit `T' `xvars' `iwgt' if `touse', base(`tbase') `psopts'
        local q
        foreach l of local tlevels {
            tempvar ql
            local q `q' `ql'
            qui predict double `ql' if `touse', pr outcome(`l')
            if "`asbalanced'"!="" {
                qui replace `ipw' = 1 / `ql' if `T'==`l' & `touse'
            }
            else {
                su `l'.`T' `awgt' if `touse', meanonly
                qui replace `ipw' = r(mean) / `ql' if `T'==`l' & `touse'
            }
        }
        if `vceadj' {
            tempname tV
            _get_V `tV'
            mata: _mktmpnames("sc", st_numscalar("e(k_eq)"))
            qui predict double `sc' if e(sample), scores
        }
    }
    // - ologit, gologit
    else if inlist("`psmethod'","ologit","gologit") {
        if "`noisily'`dots'"=="" di "..." _c
        if "`psmethod'"=="gologit" local pscmd gologit2
        else                       local pscmd `psmethod'
        `qui' `pscmd' `T' `xvars' `iwgt' if `touse', `psopts'
        if `vceadj' {
            tempvar q0 q1
            qui gen byte `q0' = 1 if `touse'
            qui gen byte `q1' = .
            local q `q0' `q1'
        }
        tempvar ql
        local j 0
        foreach l of local tlevels {
            local ++j
            qui predict double `ql' if `touse', pr outcome(`l')
            if "`asbalanced'"!="" {
                qui replace `ipw' = 1 / `ql' if `T'==`l' & `touse'
            }
            else {
                su `l'.`T' `awgt' if `touse', meanonly
                qui replace `ipw' = r(mean) / `ql' if `T'==`l' & `touse'
            }
            if `vceadj' {
                if `j'>1 {
                    qui replace `q0' = `q1' if `T'>=`l' & `touse'
                }
                if `j'<`tk' {
                    qui replace `q1' = `q0' - `ql' if `T'>=`l' & `touse'
                }
                else {
                    qui replace `q1' = 0 if `T'>=`l' & `touse'
                }
            }
            drop `ql'
        }
        if `vceadj' {
            tempname tV
            _get_V `tV'
            mata: _mktmpnames("sc", st_numscalar("e(k_eq)"))
            qui predict double `sc' if e(sample), scores
        }
    }
    // - cologit
    else if "`psmethod'"=="cologit" {
        tempvar tmpT q0 q1 ql
        qui gen byte `tmpT' = .
        qui gen byte `q0' = 1 if `touse'
        qui gen byte `q1' = .
        local q `q0' `q1'
        local tV
        local sc
        local j 0
        foreach l of local tlevels {
            if "`noisily'"=="" {
                if "`dots'"=="" di "." _c
            }
            else               di _n as res "==> level `l'"
            local ++j
            if `j'>1 {
                qui replace `q0' = `q1' if `T'>=`l' & `touse'
            }
            if `j'<`tk' {
                qui replace `tmpT' = `T'>`l' & `touse'
                `qui' logit `tmpT' `xvars' `iwgt' if `touse', `psopts'
                qui predict double `ql' if `touse', pr
                qui replace `q1' = `ql' if `T'>=`l' & `touse'
                drop `ql'
            }
            else {
                qui replace `q1' = 0 if `T'>=`l' & `touse'
            }
            qui gen double `ql' = `q0' - `q1' if `T'==`l' & `touse'
            capt assert (`ql'>=0) if `ql'<.
            if _rc==1 exit _rc
            if _rc {
                if "`noisily'`dots'"=="" di ""
                di as txt "negative probabilities encountered in level `l'; reset to 0"
                qui replace `ql' = 0 if `ql'<0 & `T'==`l' & `touse'
            }
            if "`asbalanced'"!="" {
                qui replace `ipw' = 1 / `ql' if `T'==`l' & `touse'
            }
            else {
                su `l'.`T' `awgt' if `touse', meanonly
                qui replace `ipw' = r(mean) / `ql' if `T'==`l' & `touse'
            }
            drop `ql'
            if `vceadj' & `j'<`tk' {
                tempname tVl scl
                local tV `tV' `tVl'
                local sc `sc' `scl'
                _get_V `tVl'
                qui predict double `scl' if e(sample), score
            }
        }
    }
    if "`noisily'`dots'"=="" di as txt " done)"
    
    // analyze IPWs
    capt assert (`ipw'<.) if `touse'
    if _rc==1 exit _rc
    if _rc {
        qui replace `ipw' = 1 if `ipw'>=. & `touse'
        di as err "Warning: missing values in IPW reset to 1"
    }
    qui tabstat `ipw' `awgt' if `touse', by(`T') nototal save ///
        stats(count mean sum min max cv) 
    tempname tmp ipwstats
    forv i = 1/`tk' {
        mat `tmp' = r(Stat`i')
        mat coln `tmp' = `"`r(name`i')'"'
        mat `ipwstats' = nullmat(`ipwstats') \ `tmp''
    }
    
    // estimate marginal ORs and collect results
    if "`weight'"!="" local wipw "(`wvar') * `ipw'"
    else              local wipw "`ipw'"
    if !`vceadj' {
        local vceopt `vce'
        if `"`vceopt'"'=="" local vceopt vce(robust)
    }
    else local vceopt
    logit `Y' `texpand' [iw=`wipw'] if `touse', noheader notable `vceopt' /*
        */ `constant' `mlopts'
    tempname b
    mat `b' = e(b)
    mat coleq `b' = "`depvar'"
    _ms_build_info `b' [aw=`wipw'] if `touse'
    local escalars N_cds N_cdf k df_m r2_p ll ll_0 ic rc converged
    foreach l of local escalars {
        tempname e_`l'
        scalar `e_`l'' = e(`l')
    }
    
    // VCE
    tempname V chi2 chi2_p
    mat `V' = e(V)
    mat coleq `V' = "`depvar'"
    mat roweq `V' = "`depvar'"
    if !`vceadj' {
        local vce `"`e(vce)'"'
        local vcetype `"`e(vcetype)'"'
        capt confirm matrix e(V_modelbased)
        if _rc==1 exit _rc
        if _rc==0 {
            tempname V_modelbased
            mat `V_modelbased' = e(V_modelbased)
            mat coleq `V_modelbased' = "`depvar'"
            mat roweq `V_modelbased' = "`depvar'"
        }
        else local V_modelbased `V'
        scalar `chi2' = e(chi2)
        scalar `chi2_p' = e(p)
        local chi2type `"`e(chi2type)'"'
    }
    else {
        tempname V_modelbased
        mat `V_modelbased' = `V'
    }
    if `vceadj' | "`ifgenerate'"!="" {
        mata: _mktmpnames("IFs", cols(st_matrix("e(b)")))
        _predict_logit_IF "`IFs'" "`touse'" "`texpand'" "`V_modelbased'"
        if `vceadj' {
            mata: _adjust_IFs()
            if "`clustvar'"!="" local vceopt `vce'
            else                local vceopt
            qui total `IFs' `twgt' if `touse', `vceopt'
            if "`clustvar'"!="" {
                local vce `"`e(vce)'"'
                local vcetype `"`e(vcetype)'"'
            }
            else {
                local vce robust
                local vcetype Robust
            }
            mata: st_replacematrix(st_local("V"), st_matrix("e(V)"))
        }
        else { // (unadjusted IFs)
            foreach IF of local IFs {
                qui replace `IF' = `ipw' * `IF' if `touse'
            }
        }
    }
    local rank = e(rank)
    local N_clust = e(N_clust)
    
    // post results
    eret post `b' `V' [`weight'`exp'], depname(`depvar') obs(`N') esample(`touse')
    foreach l of local escalars {
        ereturn scalar `l' = `e_`l''
    }
    eret local cmd          "ipwlogit"
    eret local title        "Marginal logistic regression"
    eret local psmethod     "`psmethod'"
    eret local psopts       `"`psopts'"'
    eret local mlopts       `"`mlopts'"'
    eret local tvar         "`tvar'"
    eret local tname        "`tname'"
    eret local ttype        "`ttype'"
    eret local tbase        "`tbase'"
    eret local indepvars    "`indepvars'"
    eret local asbaanced    "`asbalanced'"
    eret scalar k_eq        = 1
    eret scalar tk          = `tk'
    eret matrix ipwstats    = `ipwstats'
    eret local tlevels      "`tlevels'"
    if "`ttype'"!="factor" {
        if "`ttype'"!="discrete" eret scalar bins = `bins'
        eret matrix at      = `AT'
    }
    eret scalar rank        = `rank'
    eret local vce          `"`vce'"'
    eret local vcetype      `"`vcetype'"'
    eret local vceadjust    "`novceadjust'"
    eret local clustvar     "`clustvar'"
    if `"`clustvar'"'!="" {
        eret scalar N_clust = `N_clust'
    }
    capt confirm matrix `V_modelbased'
    if _rc==1 exit _rc
    else if _rc==0 {
        eret matrix V_modelbased = `V_modelbased'
    }
    if `vceadj' {
        qui test [#1]
        eret scalar chi2 = r(chi2)
        eret scalar p    = r(p)
        local chi2type "Wald"
    }
    else {
        eret scalar chi2        = `chi2'
        eret scalar p           = `chi2_p'
        eret local chi2type     `"`chi2type'"'
    }
    
    // generate
    if "`generate'"!="" {
        capt confirm variable `generate', exact
        if _rc==1 exit _rc
        if _rc==0 drop `generate'
        lab var `ipw' "Inverse probability weights"
        rename `ipw' `generate'
        eret local generate `"`generate'"'
    }
    if "`tgenerate'"!="" {
        capt confirm variable `tgenerate', exact
        if _rc==1 exit _rc
        if _rc==0 drop `tgenerate'
        lab var `T' "Binned treatment"
        rename `T' `tgenerate'
        eret local tgenerate `"`tgenerate'"'
    }
    if "`ifgenerate'"!="" {
        local coln: coln e(b)
        foreach v of local ifgenerate {
            gettoken IF IFs : IFs
            if "`IF'"=="" continue, break
            gettoken nm coln : coln
            capt confirm variable `v', exact
            if _rc==1 exit _rc
            if _rc==0 drop `v'
            lab var `IF' "IF of _b[`nm']"
            rename `IF' `v'
        }
        eret local ifgenerate "`ifgenerate'"
    }
end

program _parse_psmethod
    syntax [, Logit Mlogit Ologit GOlogit COlogit psopts(str) ]
    local psmethod `logit' `mlogit' `ologit' `gologit' `cologit'
    if `: list sizeof psmethod'>1 {
        di as err "only one method allowed in {bf:psmethod()}"
        exit 198
    }
    if "`psmethod'"=="gologit" {
        capt which gologit2
        if _rc==1 exit _rc
        if _rc {
            di as err "{bf:gologit2} is required; type {stata ssc install gologit2}"
            exit 498
        }
        _parse_psmethod_gologit_opts, `psopts'
    }
    c_local psmethod `psmethod'
end

program _parse_psmethod_gologit_opts
    syntax [, link(passthru) * ]
    if `"`link'"'!="" {
        di as err "option {bf:link()} of {bf:gologit2} not supported"
        exit 198
    }
end

program _collect_fvinfo
    local levels
    local baselevel
    foreach term of local 0 {
         _ms_parse_parts `term'
         if `"`r(type)'"'!="factor" {
             di as err "invalid specificaton of {it:tvar}"
             exit 198
         }
         if "`name'"!="" {
             if `"`r(name)'"'!="`name'" {
                 di as err "invalid specificaton of {it:tvar}"
                 exit 198
             }
         }
         local name `"`r(name)'"'
         local levels `levels' `r(level)'
         if r(omit) local baselevel `baselevel' `r(level)'
    }
    local nbase: list sizeof baselevel
    if `nbase'>1 {
        di as err "invalid specificaton of {it:tvar}"
        exit 198
    }
    if `nbase'==0 {
        // use first level as base in PS models if no base level specified
        local baselevel: word 1 of `levels'
    }
    c_local tname     "`name'"
    c_local tk        : list sizeof levels // assumed unique and sorted
    c_local tlevels   "`levels'"
    c_local tbase     "`baselevel'"
end

program _get_V
    capt confirm matrix e(V_modelbased)
    if _rc==1 exit _rc
    if _rc local Vname V
    else   local Vname V_modelbased
    mat `0' = e(`Vname')
end

program _predict_logit_IF
    args IFs touse x V
    tempvar esample sc
    qui gen byte `esample' = e(sample)==1
    qui predict double `sc' if `esample', score
    mata: _predict_logit_IF(st_local("IFs"), "`esample'", st_local("x"), ///
        "`V'", st_local("sc"))
    // reset IFs to 0 outside estimation sample (if necessary)
    capt assert (`esample') if `touse'
    if _rc==1 exit _rc
    if _rc==0 exit
    foreach IF of local IFs {
        qui replace `IF' = 0 if (!`esample') & (`touse')
    }
end

version 11
mata:
mata set matastrict on

void _mkgroups(string scalar x, string scalar g, string scalar touse,
    string scalar levels, string scalar at)
{
    real scalar    k
    real colvector X, p, l
    
    X = st_data(., x, touse)
    p = order(X, 1)
    p[p] = __mkgroups(X[p], l=.)
    k = rows(l)
    if (k>800) {
        stata(`"di as err "too many levels""')
        exit(498)
    }
    st_matrix(at, l')
    st_matrixcolstripe(at, (J(k,1,""), strofreal(1::k)))
    st_local(levels, invtokens(strofreal(1..k)))
    st_store(., st_addvar(k>100 ? "int" : "byte", g), touse, p)
}

real colvector __mkgroups(real colvector X, real matrix l)
{
    real scalar    r
    real colvector p
    
    r = rows(X)
    p = (X :!= (X[r] \ X[|1\r-1|]))
    l = select(X, p)
    return(runningsum(p))
}

void _mktmpnames(string scalar nm, real scalar k)
{
    st_local(nm, invtokens(st_tempname(k)))
}

void _predict_logit_IF(string scalar IFs, string scalar touse, string scalar xvars, 
    string scalar V, string scalar scores)
{
    real matrix sc, X, IF
    
    st_view(sc=., ., scores, touse)
    if (xvars!="") st_view(X=., ., xvars, touse)
    st_view(IF=., ., st_addvar("double", tokens(IFs)), touse)
    _fillin_IF(IF, sc, X, st_matrix(V), _predict_IF_eqinfo(V))
}

void _fillin_IF(real matrix IF, real matrix sc, real matrix X, real matrix V,
    real matrix info)
{
    real scalar j
    
    // Step 1: compute score*X
    for (j=rows(info); j; j--) {
        IF[|1,info[j,1] \ .,info[j,2]|] = _y_times_X(sc[,j], X, info[j,])
    }
    // Step 2: multiply by modelbased V
    IF[.,.] = IF * V'
}

real matrix _predict_IF_eqinfo(string scalar b)
{
    real scalar      j
    string scalar    cns
    real matrix      info
    string matrix    stripe
    
    stripe = st_matrixcolstripe(b)
    info   = panelsetup(stripe, 1)
    info   = info, J(rows(info), 1, 0)  // (from, to, hascons)
    for (j=rows(info); j; j--) {
        cns = stripe[info[j,2],2]       // last element in eq
        if (strpos(cns, ".")) cns = substr(cns, strpos(cns,".")+1, .)
        if (cns=="_cons") info[j,3] = 1 // eq has no constant
    }
    return(info)
}

real matrix _predict_IF(real matrix sc, real matrix X, real matrix V,
    real matrix info)
{
    real matrix IF

    IF = J(rows(sc), cols(V), .)
    _fillin_IF(IF, sc, X, V, info)
    return(IF)
}

real matrix _y_times_X(real colvector y, real matrix X, real rowvector info)
{
    if (info[3]) {
        if (info[1]==info[2]) return(y)  // constant only
        return((y:*X, y))                // X and constant
    }
    return(y:*X)                         // X only
}

void _adjust_IFs()
{
    real scalar      tk, tbase, l, j
    string scalar    psm, touse
    real rowvector   tl
    string rowvector sc, tV
    real colvector   T, ipw, w, p
    real matrix      q, IF, tX, tIF, tAdj, dw, info, sub
    
    // read settings and main data
    psm   = st_local("psmethod")
    tk    = strtoreal(st_local("tk"))
    tl    = strtoreal(tokens(st_local("tlevels")))
    tbase = strtoreal(st_local("tbase"))
    touse = st_local("touse")
    T     = st_data(., st_local("T"),   touse)
    ipw   = st_data(., st_local("ipw"), touse)
    w     = st_local("weight")!="" ? st_data(., st_local("wvar"), touse) : 1
    st_view(tX=., ., st_local("xvars"), touse)
    sc    = tokens(st_local("sc"))
    tV    = tokens(st_local("tV"))
    st_view(IF=., ., st_local("IFs"), touse)
    tAdj  = J(rows(IF), cols(IF), 0)
    // compute adjustment term
    // - logit
    if (psm=="logit") {
    //  + case 1: single logit
        if (tk==2 & st_local("nobinary")=="") {
            // obtain treatment model IFs
            info = _predict_IF_eqinfo(tV)
            tIF = _predict_IF(_get_sc(sc, touse), tX, st_matrix(tV), info)
            // compute adjustment
            q = st_data(., st_local("q"), touse)
            dw = w :* ((-1):^(T:!=tbase) :* (1:-q) :* ipw)
            dw = _y_times_X(dw, tX, info)
            for (j=cols(IF); j; j--) tAdj[,j] = tIF * quadcolsum(IF[,j] :* dw)'
        }
    //  + case 2: individual logit for each level
        else {                     
            q = st_data(., st_local("q"), touse)
            for (l=1; l<=tk; l++) {
                // obtain treatment model IFs
                info = _predict_IF_eqinfo(tV[l])
                tIF = _predict_IF(_get_sc(sc[l], touse), tX, st_matrix(tV[l]), info)
                // compute adjustment
                p = selectindex(T:==tl[l])
                dw = (w==1 ? w : w[p]) :* (q[p] :- 1) :* ipw[p]
                dw = _y_times_X(dw, tX[p,], info)
                for (j=cols(IF); j; j--) {
                    tAdj[,j] = tAdj[,j] + tIF * quadcolsum(IF[p,j] :* dw)'
                }
            }
        }
    }
    // - cologit
    else if (psm=="cologit") {
        q = st_data(., st_local("q"), touse)
        q = q, rowmax((J(rows(q),1,0), q[,1]-q[,2]))  // T>=t, T>t, T=t
        for (l=1; l<tk; ) {
            // obtain treatment model IFs
            info = _predict_IF_eqinfo(tV[l])
            tIF = _predict_IF(_get_sc(sc[l], touse), tX, st_matrix(tV[l]), info)
            // adjustment related to Pr(T>t)
            if (l==1) p = selectindex(T:==tl[l])
            dw = (w==1 ? w : w[p]) :* q[p,2] :* (1 :- q[p,2]) :/ q[p,3] :* ipw[p]
            dw = _y_times_X(dw, tX[p,], info)
            for (j=cols(IF); j; j--) {
                tAdj[,j] = tAdj[,j] + tIF * quadcolsum(IF[p,j] :* dw)'
            }
            // adjustment related to Pr(T>=t)
            l++
            p = selectindex(T:==tl[l])
            dw = (w==1 ? w : w[p]) :* q[p,1] :* (1 :- q[p,1]) :/ q[p,3] :* ipw[p]
            dw = _y_times_X(dw, tX[p,], info)
            for (j=cols(IF); j; j--) {
                tAdj[,j] = tAdj[,j] - tIF * quadcolsum(IF[p,j] :* dw)'
            }
        }
    }
    // - mlogit, ologit, gologit
    else {
        // obtain treatment model IFs
        info = _predict_IF_eqinfo(tV)
        tIF = _predict_IF(_get_sc(sc, touse), tX, st_matrix(tV), info)
        // compute adjustment
    // - mlogit
        if (psm=="mlogit") {
            for (l=1; l<=tk; l++) {
                q = st_data(., tokens(st_local("q"))[l], touse)
                dw = w :* ((q :- 0:^(T:!=tl[l])) :* ipw)
                dw = _y_times_X(dw, tX, info[l,])
                sub = (1, info[l,1]) \ (., info[l,2])
                for (j=cols(IF);j;j--) {
                    tAdj[,j] = tAdj[,j] + tIF[|sub|] * quadcolsum(IF[,j] :* dw)'
                }
            }
        }
    // - ologit, gologit
        else {
            q = st_data(., st_local("q"), touse)
            q = q, rowmax((J(rows(q),1,0), q[,1]-q[,2]))  // T>=t, T>t, T=t
    // - ologit
            if (psm=="ologit") {
                // main equation
                dw = w :* (q[,2] :* (1 :- q[,2]) - q[,1] :* (1 :- q[,1])) /// 
                    :/ q[,3] :* ipw
                dw = _y_times_X(dw, tX, info[1,])
                sub = (1, info[1,1]) \ (., info[1,2])
                for (j=cols(IF); j; j--) {
                    tAdj[,j] = tAdj[,j] + tIF[|sub|] * quadcolsum(IF[,j] :* dw)'
                }
                for (l=1; l<tk; ) {
                    sub = info[l+1,1]
                    // adjustment related to Tau_t
                    if (l==1) p = selectindex(T:==tl[l])
                    dw = (w==1 ? w : w[p]) :* q[p,2] :* (1 :- q[p,2]) :/ q[p,3] :* ipw[p]
                    for (j=cols(IF); j; j--) {
                        tAdj[,j] = tAdj[,j] - tIF[,sub] * quadcolsum(IF[p,j] :* dw)'
                    }
                    // adjustment related to Tau_t-1
                    l++
                    p = selectindex(T:==tl[l])
                    dw = (w==1 ? w : w[p]) :* q[p,1] :* (1 :- q[p,1]) :/ q[p,3] :* ipw[p]
                    for (j=cols(IF); j; j--) {
                        tAdj[,j] = tAdj[,j] + tIF[,sub] * quadcolsum(IF[p,j] :* dw)'
                    }
                }
            }
    // - gologit
            else if (psm=="gologit") {
                for (l=1; l<tk; ) {
                    sub = (1, info[l,1]) \ (., info[l,2])
                    // adjustment related to Pr(T>t)
                    if (l==1) p = selectindex(T:==tl[l])
                    dw = (w==1 ? w : w[p]) :* q[p,2] :* (1 :- q[p,2]) :/ q[p,3] :* ipw[p]
                    dw = _y_times_X(dw, tX[p,], info[l,])
                    for (j=cols(IF); j; j--) {
                        tAdj[,j] = tAdj[,j] + tIF[|sub|] * quadcolsum(IF[p,j] :* dw)'
                    }
                    // adjustment related to Pr(T>=t)
                    l++
                    p = selectindex(T:==tl[l])
                    dw = (w==1 ? w : w[p]) :* q[p,1] :* (1 :- q[p,1]) :/ q[p,3] :* ipw[p]
                    dw = _y_times_X(dw, tX[p,], info[l-1,])
                    for (j=cols(IF); j; j--) {
                        tAdj[,j] = tAdj[,j] - tIF[|sub|] * quadcolsum(IF[p,j] :* dw)'
                    }
                }
            }
            else exit(error(498))   // cannot be reached
        }
    }
    // adjust IFs
    IF[.,.] = ipw:*IF + tAdj
}

real matrix _get_sc(string rowvector sc, string scalar touse)
{
    return(editmissing(st_data(., sc, touse), 0))
}

end

