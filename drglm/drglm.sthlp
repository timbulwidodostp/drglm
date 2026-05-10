{smcl}
{* *! version 1.0.0 25sep12}{...}
{cmd:help drglm}{right: ({browse "http://www.stata-journal.com/article.html?article=st0290":SJ13-1: st0290})}
{hline}

{title:Title}

{p2colset 5 13 10 5}{...}
{p2col :{hi:drglm} {hline 1}}Doubly robust estimation in generalized linear models{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 10 2}
{cmd:drglm}
{depvar}
{it:expvar}
{ifin}
[{cmd:,} {opt m:ain(varlist)} 
{opt o:utcome(varlist)} 
{opt e:xposure(varlist)} 
{opt ol:ink(linkname)}
{opt el:ink(linkname)}
{opt l:evel(#)}
{opt obe}
{opt ebe}
{opt ef:orm}
{opth vce(vcetype)}]

{pstd}{cmd:by} and {cmd:mi estimate} are allowed; see {help prefix}.


{title:Description}

{pstd}{cmd:drglm} provides doubly robust estimation for the main model
in generalized linear models.  The {it:expvar} (exposure, treatment,
predictor, or covariate) must be numerical.

{pstd}After {cmd:drglm} estimation, one can use postestimation commands
such as {helpb test}, {helpb testparm}, {helpb lincom}, and 
{helpb predictnl}.


{title:Options}

{phang}{opt main(varlist)} determines which variables are used in the
main model part of the estimator.  The constant 1 is always added to
{opt main(varlist)}.  Then each variable in {opt main(varlist)} is
multiplied by {it:expvar} and saved in the current dataset.

{phang}{opt outcome(varlist)} determines which variables are used in the
outcome model part of the estimator.  The constant 1 is always added to
{opt outcome(varlist)}.

{phang}{opt exposure(varlist)} determines which variables are used in
the exposure model part of the estimator.  The constant 1 is always
added to {opt exposure(varlist)}.

{phang}{opt olink(linkname)} specifies the link function of the outcome
model ({cmd:identity}, {cmd:logit}, {cmd:log}).  The default is
{cmd:olink(identity)}.  See table below.  If {cmd:olink(logit)} is
specified, {it:expvar} can take on only two values (either 0 or 1).

{phang}{opt elink(linkname)} specifies the link function of the exposure
model ({cmd:identity}, {cmd:logit}, {cmd:log}).  The default is
{cmd:elink(identity)}.  See table below.

{phang}{opt level(#)} specifies the confidence level, as a percentage,
for confidence intervals.  The default is {cmd:level(95)} or as set by
{helpb set level}.

{phang}{opt obe} specifies the outcome-based estimation.

{phang}{opt ebe} specifies the exposure-based estimation.

{phang}{opt eform} reports coefficient estimates as {cmd:exp(b)} rather
than as {cmd:b}.

{phang}{opt vce(vcetype)} specifies the type of standard error reported.
{it:vcetype} may be {cmd:robust}, {opt cl:uster} {it:clustvar}, 
{opt boot:strap}, or {opt jack:knife}.  The default is
{cmd:vce(robust)}.


{title:Possible combinations of link functions}

        {bf:olink()           elink()}
	{hline 30}
	{cmd:identity} 		{cmd:identity}
	{cmd:identity} 		{cmd:log}
	{cmd:identity} 		{cmd:logit}
	{cmd:log} 		{cmd:identity}
	{cmd:log} 		{cmd:log}
	{cmd:log} 		{cmd:logit}
	{cmd:logit} 		{cmd:logit}


{title:Examples}

{phang2}{cmd:.} {bf:{stata "webuse lbw"}}{p_end}

{pstd}Doubly robust estimation for continuous outcomes{p_end}
{phang2}{cmd:.} {bf:{stata "xi: drglm  bwt smoke, outcome(age lwt i.race)"}}{p_end}

{pstd}Doubly robust estimation for binary outcomes{p_end}
{phang2}{cmd:.} {bf:{stata "xi: drglm low smoke, outcome(age lwt i.race) olink(logit) eform"}}{p_end}


{title:Saved results}

{pstd}
{cmd:drglm} saves the following in {cmd:e()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(rank)}}rank of {cmd:e(V)}{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:drglm}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(vcetype)}}title used to label Std. Err.{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(olink)}}link function of the outcome model{p_end}
{synopt:{cmd:e(elink)}}link function of the exposure model{p_end}
{synopt:{cmd:e(estimator)}}type of estimator ({cmd:dr}, {cmd:obe}, or {cmd:ebe}){p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}


{title:Authors}

{pstd}Nicola Orsini{p_end}
{pstd}Unit of Biostatistics and Unit of Nutritional Epidemiology{p_end}
{pstd}Institute of Environmental Medicine{p_end}
{pstd}Karolinska Institutet{p_end}
{pstd}Stockholm, Sweden{p_end}
{pstd}nicola.orsini@ki.se{p_end}

{pstd}Rino Bellocco{p_end}
{pstd}Department of Medical Epidemiology and Biostatistics{p_end}
{pstd}Karolinska Institutet{p_end}
{pstd}Stockholm, Sweden{p_end}
{pstd}rino.bellocco@ki.se{p_end}

{pstd}Arvid Sj{c o:}lander{p_end}
{pstd}Department of Medical Epidemiology and Biostatistics{p_end}
{pstd}Karolinska Institutet{p_end}
{pstd}Stockholm, Sweden{p_end}
{pstd}arvid.sjolander@ki.se{p_end}
	

{title:Also see}

{p 4 14 2}Article:  {it:Stata Journal}, volume 13, number 1: {browse "http://www.stata-journal.com/article.html?article=st0290":st0290}

{p 5 14 2}Manual:  {manlink R glm}, {manlink R gmm}
{p_end}
