
program define gaci_core, rclass
    version 17
    syntax , OUTcome(varname numeric) EXPosure(varname numeric) Z(varlist numeric) ///
        [ U(varlist numeric) SCALE(string) ANGTHResh(real 15) KAPPA0(real 0.05) ]

    if "`scale'"=="" local scale "logit"
    tempvar touse
    marksample `touse'
    quietly keep if `touse'

    capture drop _gaciX
    quietly summarize `exposure'
    gen double _gaciX = (`exposure' - r(mean)) / r(sd)

    local Zc
    local zi = 0
    foreach v of local z {
        local ++zi
        capture drop _gaciZ`zi'
        quietly summarize `v'
        gen double _gaciZ`zi' = (`v' - r(mean)) / r(sd)
        local Zc `Zc' _gaciZ`zi'
    }
    foreach v of local Zc {
        capture drop _gaciXZ_`v'
        gen double _gaciXZ_`v' = _gaciX*`v'
    }
    local XZlist
    unab XZlist : _gaciXZ_*

    if "`scale'"=="logit" {
        glm `outcome' _gaciX `Zc' `XZlist', family(binomial) link(logit) vce(robust)
    }
    else if "`scale'"=="log" {
        glm `outcome' _gaciX `Zc' `XZlist', family(binomial) link(log) vce(robust)
    }
    else if "`scale'"=="identity" {
        glm `outcome' _gaciX `Zc' `XZlist', family(gaussian) link(identity) vce(robust)
    }
    else {
        di as err "scale() must be logit | log | identity"
        exit 198
    }

    matrix V = e(V)
    matrix I = invsym(V)
    scalar det_g = det(I)
    mata: eig = eigenvalue(st_matrix("I"))
    mata: st_numscalar("min_eig", min(eig))
    scalar metric_ok = (det_g>1e-8 & min_eig>1e-8)

    mata:
        S = GACI()
        S.b     = st_matrix("e(b)")
        S.names = st_matrixcolstripe("e(b)")
        S.link  = st_local("scale")
        S.xname = "_gaciX"
        string scalar zlist = st_local("Zc")
        string rowvector zt = tokens(zlist)
        S.znames = zt
    end

    local k = wordcount("`Zc'")
    mata: p0 = (0, J(1, `k', 0))
    mata: plo = (0, J(1, `k', -1))
    mata: phi = (0, J(1, `k', 1))

    mata: kappa = gaci_kappa(S, p0, 1e-4)
    scalar kappa = kappa

    mata:
        real colvector v_lo = gaci_transport_hess(S, plo, plo, 1e-4)
        real colvector v_hi = gaci_transport_hess(S, phi, phi, 1e-4)
        real colvector v_tr = gaci_transport_hess(S, plo, phi, 1e-4)
        I = st_matrix("I")
        theta = gaci_angle(v_tr, v_hi, I)
        real colvector v_back = gaci_transport_hess(S, phi, plo, 1e-4)
        hol = gaci_angle(v_lo, v_back, I)
    end
    scalar theta    = theta
    scalar holonomy = hol

    quietly regress `outcome' _gaciX `Zc'
    scalar beta_adjZ = _b[_gaciX]
    mata: grad_lo = gaci_grad(S, plo, 1e-4)
    mata: st_numscalar("beta_causal", grad_lo[1,1])
    scalar dprojZ = 1 - (beta_adjZ*beta_causal)/(sqrt(beta_adjZ^2)*sqrt(beta_causal^2))

    scalar dprojZU = .
    if "`u'" != "" {
        quietly regress `outcome' _gaciX `Zc' `u'
        if !_rc {
            scalar beta_adjZU = _b[_gaciX]
            scalar dprojZU = 1 - (beta_adjZU*beta_causal)/(sqrt(beta_adjZU^2)*sqrt(beta_causal^2))
        }
    }

    return scalar kappa     = kappa
    return scalar theta     = theta
    return scalar dprojZ    = dprojZ
    return scalar dprojZU   = dprojZU
    return scalar holonomy  = holonomy
    return scalar metric_ok = metric_ok
    return scalar det_g     = det_g
    return scalar min_eig   = min_eig
end

program define gaci_run, rclass
    version 17
    syntax , NAME(string) OUTcome(varname numeric) EXPosure(varname numeric) Z(varlist numeric) ///
        [ U(varlist numeric) SCALE(string) ANGTHResh(real 15) KAPPA0(real 0.05) SAVE(string "gaci_results.dta") BSREPS(integer 0) ]

    quietly gaci_core, outcome(`outcome') exposure(`exposure') z(`z') u(`u') scale(`scale') angthresh(`angthresh') kappa0(`kappa0')
    foreach s in kappa theta dprojZ dprojZU holonomy metric_ok det_g min_eig {
        scalar `s' = r(`s')
    }

    tempname bV
    if `bsreps' > 0 {
        noisily di as text "Bootstrapping (" `bsreps' " reps) ..."
        bootstrap kappa=r(kappa) theta=r(theta), reps(`bsreps') nowarn nodots: ///
            gaci_core, outcome(`outcome') exposure(`exposure') z(`z') u(`u') scale(`scale') angthresh(`angthresh') kappa0(`kappa0')
        capture scalar kappa_se = _se[kappa]
        capture scalar theta_se = _se[theta]
    }

    tempname mem
    capture confirm file "`save'"
    if _rc {
        postfile `mem' str32 name strL outcome strL exposure strL z strL u str8 scale ///
            double kappa theta dprojZ dprojZU holonomy byte metric_ok double det_g min_eig long N ///
            double kappa_se theta_se using "`save'", replace
    }
    else {
        postfile `mem' str32 name strL outcome strL exposure strL z strL u str8 scale ///
            double kappa theta dprojZ dprojZU holonomy byte metric_ok double det_g min_eig long N ///
            double kappa_se theta_se using "`save'", append
    }
    quietly count
    local N = r(N)
    post `mem' ("`name'") ("`outcome'") ("`exposure'") ("`z'") ("`u'") ("`scale'") ///
        (kappa) (theta) (dprojZ) (dprojZU) (holonomy) (metric_ok) (det_g) (min_eig) (`N') ///
        (cond(`bsreps'>0, kappa_se, .)) (cond(`bsreps'>0, theta_se, .))
    postclose `mem'

    di as text _n "Cycle: " as result "`name'" _n "  κ = " %6.4f kappa "   θ = " %6.2f theta "°   Δ_proj(Z) = " %6.3f dprojZ "   Δ_proj(Z+U) = " %6.3f dprojZU "   holonomy = " %6.2f holonomy
    if metric_ok==0 di as err "  [Warning] Fisher metric unstable: det(g)=" %9.2e det_g ", λ_min=" %9.2e min_eig
    if `bsreps' > 0 di as text "  [Bootstrap] se(κ)=" %6.4f kappa_se "  se(θ)=" %6.4f theta_se
end

// --- Mata utilities injected below ---

// ---- GACI Mata utilities (shared) ----
mata:
string scalar _b_get_idx(string scalar coef, string matrix names) {
    real scalar idx = 0
    for (i=1; i<=rows(names); i++) {
        if (names[i,2] == coef) { idx = i; break }
    }
    return(strofreal(idx))
}

real scalar _b_get(string scalar coef, real rowvector b, string matrix names) {
    real scalar idx = real(_b_get_idx(coef, names))
    if (idx==0) return(0)
    else return(b[1,idx])
}

struct GACI {
    real rowvector b
    string matrix names
    string scalar link
    string scalar xname
    string rowvector znames
};

real scalar gaci_gmu(struct GACI scalar S, real rowvector p)
{
    real scalar x = p[1,1]
    real scalar k = cols(p)-1
    real rowvector z = (k>0 ? p[1,2..cols(p)] : J(1,0,.))

    real scalar eta = _b_get("_cons", S.b, S.names) + _b_get(S.xname, S.b, S.names) * x
    for (j=1; j<=k; j++) {
        string scalar zj = S.znames[1,j]
        eta = eta + _b_get(zj, S.b, S.names) * z[1,j]
        string scalar ij = "c." + S.xname + "#c." + zj
        eta = eta + _b_get(ij, S.b, S.names) * x * z[1,j]
    }
    return(eta)
}

real rowvector gaci_grad(struct GACI scalar S, real rowvector p, real scalar h)
{
    real rowvector g = J(1, cols(p), .)
    for (i=1; i<=cols(p); i++) {
        real rowvector ei = J(1, cols(p), 0); ei[1,i]=h
        g[1,i] = (gaci_gmu(S, p:+ei) - gaci_gmu(S, p:-ei)) / (2*h)
    }
    return(g)
}

real matrix gaci_hessian(struct GACI scalar S, real rowvector p, real scalar h)
{
    real scalar n = cols(p)
    real matrix H = J(n,n,0)
    for (i=1; i<=n; i++) {
        for (j=i; j<=n; j++) {
            real rowvector ei = J(1,n,0); ei[1,i]=h
            real rowvector ej = J(1,n,0); ej[1,j]=h
            real scalar fpp = gaci_gmu(S, p:+ei:+ej)
            real scalar fpm = gaci_gmu(S, p:+ei:-ej)
            real scalar fmp = gaci_gmu(S, p:-ei:+ej)
            real scalar fmm = gaci_gmu(S, p:-ei:-ej)
            H[i,j] = (fpp - fpm - fmp + fmm)/(4*h*h)
            H[j,i] = H[i,j]
        }
    }
    return(H)
}

real scalar gaci_kappa(struct GACI scalar S, real rowvector p, real scalar h)
{
    real matrix H = gaci_hessian(S, p, h)
    real rowvector grad = gaci_grad(S, p, h)
    real rowvector gz = grad; gz[1,1]=0
    real scalar ngz = norm(gz)
    if (ngz>0) gz = gz/ngz
    real rowvector ex = J(1, cols(p), 0); ex[1,1]=1
    real matrix B = (ex \ gz)
    real matrix Hpi = B * H * B'
    return( sqrt(sum(Hpi:^2)) )
}

real scalar gaci_angle(real colvector v, real colvector w, real matrix S)
{
    real scalar num = v' * S * w
    real scalar dv = sqrt(v' * S * v)
    real scalar dw = sqrt(w' * S * w)
    if (dv==0 | dw==0) return(.)
    real scalar c = num / (dv * dw)
    if (c>1) c=1
    if (c<-1) c=-1
    return(acos(c)*180/pi())
}

real colvector gaci_transport_hess(struct GACI scalar S, real rowvector p, real rowvector q, real scalar h)
{
    real rowvector gp = gaci_grad(S, p, h)
    real rowvector gq = gaci_grad(S, q, h)
    real matrix Hp = gaci_hessian(S, p, h)
    real matrix Hq = gaci_hessian(S, q, h)
    real rowvector d = q - p
    real rowvector Havg = 0.5*(Hp:+Hq)
    real rowvector upd = d * Havg'
    real scalar v_x = gp[1,1] + upd[1,1]
    real scalar v_zmag = (norm(gp[1,2..cols(gp)]) + norm(gq[1,2..cols(gq)]))/2
    real colvector v2 = (v_x \ v_zmag)
    return(v2)
}
end

