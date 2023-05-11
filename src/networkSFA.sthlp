{smcl}
{vieweralsosee "[R] help" "help help"}{...}
{viewerjumpto "Syntax" "networkSFA##syntax"}{...}
{viewerjumpto "Description" "networkSFA##description"}{...}
{viewerjumpto "Options" "networkSFA##options"}{...}
{viewerjumpto "Remarks" "networkSFA##remarks"}{...}
{viewerjumpto "Examples" "networkSFA##examples"}{...}
{title:Title}

{phang}
{bf:networkSFA} {hline 2} Network two-stage stochastic frontier model with no shared inputs (i.e. All inputs are used in the first production stage)


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:networkSFA}
{case}
({depvar1} {varlist1})
({depvar2} {varlist2})
{...}
({depvarN} {varlistN})
[{list of IOs if case is IO}]
[{if}]


{title:Description}

{pstd}
{cmd:networkSFA} estimates of the network two-stage SFA model where the first production stage produces intermediate outputs which are used to produce final outputs in the second stage. There is no other inputs (rather than intermediate outputs) involved in the second stage. The command allows for the models with different data aggregations specified in our paper. The model is estimated by a multi-step procedure. In the first step, SUR or 3SLS is used to estimate slope parameters. The second step uses the method of moments (MM) to estimate distributional parameters and correct intercepts. The last step is to gauge efficiency. 

{title:Authors}

{pstd}Vo Huyen Trang Tran{p_end}
{pstd}Institute for Transport Studies, University of Leeds{p_end}
{pstd}Leeds, UK{p_end}
{pstd}v.h.tran@leeds.ac.uk{p_end}

{title:Also see}

{psee}
{space 2}Help:  {help network_postestimation}
{p_end}
