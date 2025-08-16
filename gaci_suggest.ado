
program define gaci_suggest, rclass
    version 17
    syntax [, SHOWall]
    tempname touse
    marksample `touse'
    ds, has(type numeric)
    local nvars `r(varlist)'
    local binlist
    foreach v of local nvars {
        quietly summarize `v' if `touse'
        if r(N)>0 & r(min)==0 & r(max)==1 {
            local binlist `binlist' `v'
        }
    }
    local contlist : list nvars - binlist
    di as text "---- Suggested variables ----"
    di as text "Outcome (binary Y) candidates:" _n as result "`binlist'"
    di as text "Exposure (X) candidates:" _n as result "`binlist' `contlist'"
    di as text "Covariates (Z) candidates:" _n as result "`contlist'"
    return local Ycandidates "`binlist'"
    return local Xcandidates "`binlist' `contlist'"
    return local Zcandidates "`contlist'"
end
