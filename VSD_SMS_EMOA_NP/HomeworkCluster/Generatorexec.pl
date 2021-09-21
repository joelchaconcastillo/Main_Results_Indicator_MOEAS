#Con este script se genera la lista de ejecuciones....
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my $file = "ExecutionFileDiversity_d4";
my $fout;
open($fout, '>' ,$file);
my $Path =  `cd ..; pwd`;#"/home/joel.chacon/Current/MyResearchTopics/MOEA-D-Diversity/MOEAD-DE/vsd-moead-opt";
chomp $Path;

####Realizar la búsqueda del parámetro D inicial que proporcione mejores resultados
my $PathAlgorithm = $Path;

#for(my $Df=0.1; $Df <= 0.9; $Df+=0.2)
#{
#my @Instance = ("DTLZ1", "DTLZ2", "DTLZ3", "DTLZ4", "DTLZ5", "DTLZ6", "DTLZ7");
my @DI = ("0.4");#, "0.2", "0.4", "0.6", "0.8", "1.0");
#foreach my $CR("0.2", "0.4", "0.6", "0.8", "1.0")
foreach my $CR("0.2") #, "0.6", "0.8", "1.0")
{
foreach my $Di(@DI)
{
	#foreach my $Df(("0.1","0.2","0.3","0.4")){
	foreach my $Df(("0.5")){
	#   my $Di = 0.4;
	#   my $Df = 0.5;
   my $pops=100;
   my $max_nfes=25000000;
   #   my $CR = 0.0;
   my $F = 0.75;
   my $nOffspring = 1; 
   
   my @Instance = ("DTLZ1", "DTLZ2", "DTLZ3", "DTLZ4", "DTLZ5", "DTLZ6", "DTLZ7");
            foreach(@Instance)
            {
            	my $nvar;
            	
            	for(my $nobj = 2; $nobj <=3; $nobj++)
            	{
               	   if($_ eq "DTLZ1")
            	   {
            	      $nvar=5+$nobj-1;
            	   }
            	   elsif($_ eq "DTLZ7")
            	   {
            	      $nvar=20+$nobj-1;
            	   }
            	   else
            	   {
            	      $nvar=10+$nobj-1;
            	   }
            	
            	   for(my $Sed = 1; $Sed <=35; $Sed++) ##Realizar 35 ejecuciones con distintas semilla de cada instancia..
            	   {
            	   	print $fout "~$PathAlgorithm/Ejecutable $PathAlgorithm $_ $Sed $nobj $pops $nOffspring $max_nfes $CR $F $nvar $Di $Df\n";
            	   }
            	}
            }
      @Instance = ("WFG1", "WFG2", "WFG3", "WFG4", "WFG5", "WFG6", "WFG7", "WFG8", "WFG9");
      #   @Instance = ("WFG8");
            foreach(@Instance)
            {
            	my $nvar;
            	
            	for(my $nobj = 2; $nobj <=3; $nobj++)
            	{
            	   my $k =4;# 2*($nobj-1);
            	   my $l =20;# 24-$k;
            	   $nvar=$l+$k;
    	
            	   for(my $Sed = 1; $Sed <=35; $Sed++) ##Realizar 35 ejecuciones con distintas semilla de cada instancia..
            	   {
            	   	print $fout "~$PathAlgorithm/Ejecutable $PathAlgorithm $_ $Sed $nobj $pops $nOffspring $max_nfes $CR $F $nvar $Di $Df\n";
            	   }
            	}
            }
	       	 	    @Instance = ("UF1", "UF2", "UF3", "UF4", "UF5", "UF6", "UF7");
	      	           my $nobj=2;
	       	 	       foreach(@Instance)
	       	 	       {
	       	 	       	   my $nvar=30;
	       	 	       	   for(my $Sed = 1; $Sed <=35; $Sed++) ##Realizar 35 ejecuciones con distintas semilla de cada instancia..
	       	 	       	   {
	       	 	       	   	print $fout "~$PathAlgorithm/Ejecutable $PathAlgorithm $_ $Sed $nobj $pops $nOffspring $max_nfes $CR $F $nvar $Di $Df\n";
	       	 	       	   }
	       	 	       }
   	    @Instance = ("UF8", "UF9", "UF10");
   	    $nobj=3;
   	       foreach(@Instance)
   	       {
   	       	   my $nvar=30;
   
   	       	   for(my $Sed = 1; $Sed <=35; $Sed++) ##Realizar 35 ejecuciones con distintas semilla de cada instancia..
   	       	   {
   	       	   	print $fout "~$PathAlgorithm/Ejecutable $PathAlgorithm $_ $Sed $nobj $pops $nOffspring $max_nfes $CR $F $nvar $Di $Df\n";
   	       	   }
   	       }

	       	 	    @Instance = ("IMB1", "IMB2", "IMB3", "IMB7", "IMB8", "IMB9");
	       	           $nobj=2;
	       	 	       foreach(@Instance)
	       	 	       {
	       	 	       	   my $nvar=10;
	       	 	       	   for(my $Sed = 1; $Sed <=35; $Sed++) ##Realizar 35 ejecuciones con distintas semilla de cada instancia..
	       	 	       	   {
	       	 	       	   	print $fout "~$PathAlgorithm/Ejecutable $PathAlgorithm $_ $Sed $nobj $pops $nOffspring $max_nfes $CR $F $nvar $Di $Df\n";
	       	 	       	   }
	       	 	       }
	       	    
	    @Instance = ("IMB4", "IMB5", "IMB6", "IMB10");
	    $nobj=3;
	       foreach(@Instance)
	       {
	       	   my $nvar=10;

	       	   for(my $Sed = 1; $Sed <=35; $Sed++) ##Realizar 35 ejecuciones con distintas semilla de cada instancia..
	       	   {
	       	   	print $fout "~$PathAlgorithm/Ejecutable $PathAlgorithm $_ $Sed $nobj $pops $nOffspring $max_nfes $CR $F $nvar $Di $Df\n";
	       	   }
	       }

       }
}
}
