{smcl}
{cmd:help networkSFA postestimation} {...}
{right:also see:  {help networkSFA}}
{hline}

{title:Title}

{p2colset 5 32 38 2}{...}
{p2col :{hi:networkSFA postestimation} {hline 2}}Postestimation tools for
networkSFA{p_end}
{p2colreset}{...}

{title:Description}

{pstd}
The following postestimation commands are available for {opt networkSFA}:

{synoptset 11}{...}
{p2coldent :command}description{p_end}
{synoptline}
{synopt :{helpb rfrontier postestimation##predict:predict}}predictions, residuals,
efficiency scores{p_end}
{synoptline}
{p2colreset}{...}

{marker predict}{...}
{title:Syntax for predict}

{p 8 16 2}{cmd:predict} {newvar} {help if: if} [{cmd:,} {it:statistic} eq(eqno)]

{synoptset 15 tabbed}{...}
{synopthdr :statistic}
{synoptline}
{syntab :Main}
{synopt :{opt xb}}linear prediction; the default{p_end}
{synopt :{opt residuals}}residuals{p_end}
{synopt :{opt u}}estimates of (technical or cost) inefficiency via {it:E}(u|e)
(Jondrow et al., 1982){p_end}
{synopt :{opt m}}estimates of (technical or cost) inefficiency via {it:M}(u|e)
(Jondrow et al., 1982){p_end}
{synopt :{opt jlms}}estimates of (technical or cost) efficiency via exp[-{it:E}(u|e)]
(Jondrow et al., 1982){p_end}
{synopt :{opt bc}}estimates of (technical or cost) inefficiency {it:E}[exp(-u|e)]
(Battese and Coelli, 1988){p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
eq(eqno) specifies the equation you would like to calculate a statistic.  eq(1) means the first equation of the system (the default), eq(2) means the second equation and so on.  eq(0) is only used for (in)efficiency statistics (u, m, jlms, bc). It is to specify that you would like to calculate overall (in)efficiency.

{title:Authors}

{pstd}Vo Huyen Trang Tran{p_end}
{pstd}Institute for Transport Studies, University of Leeds{p_end}
{pstd}Leeds, UK{p_end}
{pstd}v.t.tran@leeds.ac.uk{p_end}

