
program define gaci_compare
    version 17
    syntax [, USING(string) SHOW ANGTHResh(real 15) AUTOthresh ]
    if "`using'"=="" local using "gaci_results.dta"
    use "`using'", clear
    order name scale kappa theta dprojZ dprojZU holonomy det_g min_eig N kappa_se theta_se, first
    format kappa %6.4f theta %6.2f dprojZ %6.3f dprojZU %6.3f holonomy %6.2f det_g %9.2e min_eig %9.2e
    list, abbrev(20) noobs sepby(scale)

    local kappa0 = .
    if "`autothresh'" != "" {
        quietly summarize kappa, detail
        local kappa0 = r(p95)
        di as text "Auto-thresholds: κ{c 0xb0} (p95) = " %6.4f `kappa0' " ; θ{c 0xb0} = " %5.1f `angthresh' "°"
    }

    if "`autothresh'" != "" & `kappa0' < . {
        twoway ///
            (scatter theta kappa, mlabel(name) mlabpos(0)) ///
            (function y=`angthresh', range(min(kappa) max(kappa)) lpattern(dash)) ///
            (function y=x*0+., range(`kappa0' `kappa0') vert lpattern(dash)), ///
            xtitle("{&kappa} (curvature)") ytitle("{&theta} (degrees)") name(GACI_theta_kappa, replace)
    }
    else {
        twoway scatter theta kappa, mlabel(name) mlabpos(0) xtitle("{&kappa} (curvature)") ytitle("{&theta} (degrees)") name(GACI_theta_kappa, replace)
    }
end
