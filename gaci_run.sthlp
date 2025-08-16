
{smcl}
{title:Title}
{phang}
{bf:gaci_run} — Run one causal-geometry diagnostic cycle (with optional bootstrap)

{title:Syntax}
{p 8 12 2}
{cmd:gaci_run}, {opt name(string)} {opt out:come(varname)} {opt exp:osure(varname)} {opt z(varlist)} [{opt u(varlist)} {opt scale(string)} {opt angthr:esh(#)} {opt kappa0(#)} {opt save(string)} {opt bsreps(#)}]

{title:Description}
{pstd}
Fits {it:μ(x,z)} on the chosen scale (logit/log/identity), computes {it:κ = ||H_Π||_F}, transports effect vectors via a
Hessian-propagation step, and reports {it:θ}, {it:Δ} (Z and Z+U), holonomy, and metric diagnostics. Optionally bootstraps
{it:κ} and {it:θ}.

{title:Options}
{synoptset 22 tabbed}{...}
{synopt:{opt bsreps(#)}}If >0, perform bootstrap on κ and θ with {it:#} replications; returns SEs when available{p_end}
{synopt:{opt scale(string)}}{bf:logit} (default), {bf:log}, or {bf:identity}{p_end}
{synopt:{opt u(varlist)}}Optional U-proxies/negative controls for Δ_proj(Z+U){p_end}
{synopt:{opt save(string)}}Output results dataset; default {it:gaci_results.dta}{p_end}

{title:Returned scalars}
{synoptset 25}{...}
{synopt:{cmd:r(kappa) r(theta) r(dprojZ) r(dprojZU) r(holonomy)}}key diagnostics{p_end}
{synopt:{cmd:r(metric_ok) r(det_g) r(min_eig)}}metric diagnostics{p_end}

{title:Notes}
{pstd}
Hessian-propagation is consistent with the gradient field of {it:g(μ)}; thresholds should be dataset-calibrated (bootstrap).

{title:Examples}
{phang2}{cmd:. gaci_run, name(A) outcome(y) exposure(treat) z(age bmi smoke) scale(logit) bsreps(200)}
