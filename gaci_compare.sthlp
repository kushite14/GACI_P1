
{smcl}
{title:Title}
{phang}
{bf:gaci_compare} — Compare cycles and draw θ–κ with optional auto-thresholds

{title:Syntax}
{p 8 12 2}
{cmd:gaci_compare} [{opt using(filename)}] [{opt angthr:esh(#)} {opt autothresh} {opt show}]

{title:Description}
{pstd}
Displays a comparison table and generates a θ–κ scatter. With {opt autothresh}, overlays κ₀ at the 95th percentile of κ and θ₀ at {opt angthresh}.

{title:Options}
{synoptset 22}{...}
{synopt:{opt autothresh}}Auto-set κ₀ = p95(κ); draw vertical line and horizontal θ₀{p_end}
{synopt:{opt angthr:esh(#)}}θ₀ guide (degrees); default 15{p_end}
{synopt:{opt show}}Display full diagnostics table with metric details{p_end}

{title:Option Details}
{p2colset 9 28 32}{...}
{p2col:{opt show}}When specified, displays the complete comparison table including:
{p_end}
{p2col:9 28 32:}• {it:holonomy} (parallel transport inconsistency){p_end}
{p2col:9 28 32:}• {it:det_g} (determinant of Fisher metric){p_end}
{p2col:9 28 32:}• {it:min_eig} (smallest eigenvalue){p_end}
{p2col:9 28 32:}• {it:metric_ok} flag (1=stable metric){p_end}
{p2col:9 28 32:}• Bootstrap SEs (if available){p_end}
