
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="EN">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="30%" cellspacing="0" cellpadding="0">
<tr><td><a href="http://r-forge.r-project.org/"><img src="rforgelogo.png" border="0" alt="R-Forge Logo" /> </a> </td>
    <td><a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><img src="logo.png" align="left" border="0" alt="DDEPN Logo" /> </a></td></tr>
</table>

<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<!-- h1 align="left"> Welcome to the ddepn project.</h1 -->



<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description --> 
<p> This project contains the
method 'Dynamic Deterministic Effects Propagation Networks',
implemented in the R-package 'ddepn'. The package provides a network
inference engine for high-throughput genomic or proteomic expression
data, generated after external intervention (inhibitions or
stimulations). Two main parts are implemented: an MCMC network
structure sampling approach and a Genetic Algorithm network
optimisation. Further biological prior knowledge can be included using
different prior models.</p><br/>

<!-- dl: definition list, dt: definition title, dd: definition
description -->

The method is described in the following publications:<br/>
<ul>
  <li><dl>
  <dt>Article</dt>
    <dd>Bender, C., Henjes, F., Fr&ouml;hlich, H., Wiemann, S., Korf, U. &amp; Bei&szlig;barth, T.</dd>
    <dd><a href="http://www.ncbi.nlm.nih.gov/pubmed/20823327"
	   target="_blank">Dynamic deterministic effects propagation
	     networks: learning signalling pathways from longitudinal
	     protein array data.</a></dd> 
    <dd><em>Bioinformatics, </em><strong>2010</strong>, Vol. 26(18), pp. i596-i602</dd>
  </dl></li>
  <li><dl>
  <dt>PhD Thesis</dt>
    <dd>Bender, C.</dd>
    <dd><a href="http://www.dkfz.de/mga2/people/bender/files/BenderThesis.pdf" target="_blank">
    	   Systematic analysis of time resolved high-throughput data using stochastic network inference methods</a></dd> 
    <dd><em>University of Heidelberg, Combined Faculties for the Natural Sciences and for Mathematics, </em><strong>2011</strong></dd>
  </dl></li>	       
  <li><dl>	       
    <dt>Article</dt>
	  <dd>Bender, C., Heyde, S. v.d., Henjes, F., Wiemann, S.,
	  Korf, U. &amp; Bei&szlig;barth, T.</dd>
	  <dd><a href="http://www.ncbi.nlm.nih.gov/pubmed/21771315" target="_blank">Inferring signalling networks from longitudinal data
	  using sampling based approaches in the R package 'ddepn'</a>
	  </dd> 
	  <dd><em>BMC Bioinformatics, </em><b>2011,</b><i> 12:291</i></dd>
  </dl></li>
</ul>

<p> The <strong>project summary page</strong>
you can find <a href="http://<?php echo $domain; ?>/projects/<?php
echo $group_name; ?>/"><strong>here</strong></a>. Also visit the DKFZ homepage for more information on the <a href="http://www.dkfz.de/mga2/people/bender">author</a>.</p>

</body>
</html>
