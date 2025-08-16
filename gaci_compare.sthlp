
{smcl}
{title:Title}
{phang}
{bf:gaci_compare} — Compare cycles and draw θ–κ with optional auto-thresholds

{title:Syntax}
{p 8 12 2}
{cmd:gaci_compare} [{opt using(filename)}] [{opt angthr:esh(#)} {opt autothresh}]

{title:Description}
{pstd}
Displays a comparison table and generates a θ–κ scatter. With {opt autothresh}, overlays κ₀ at the 95th percentile of κ and θ₀ at {opt angthresh}.

{title:Options}
{synoptset 22}{...}
{synopt:{opt autothresh}}Auto-set κ₀ = p95(κ); draw vertical line and horizontal θ₀{p_end}
{synopt:{opt angthr:esh(#)}}θ₀ guide (degrees); default 15{p_end}
