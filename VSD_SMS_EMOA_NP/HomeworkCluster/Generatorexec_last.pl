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
my $Di = 0.4;
my $Df = 0.5;
my $pops=100;
my $max_nfes=250000*$pops;
my $CR = 0.0;
my $F = 0.75;
my @vars = ("50", "100", "250", "500");

foreach my $nvar(@vars)
{
      my @Instance = ("DTLZ1", "DTLZ2", "DTLZ3", "DTLZ4", "DTLZ5", "DTLZ6", "DTLZ7");
         foreach(@Instance)
         {
         	for(my $nobj = 2; $nobj <=3; $nobj++)
         	{
         	   for(my $Sed = 1; $Sed <=35; $Sed++) ##Realizar 35 ejecuciones con distintas semilla de cada instancia..
         	   {
         	   	print $fout "~$PathAlgorithm/Ejecutable $PathAlgorithm $_ $Sed $nobj $pops $max_nfes $CR $F $nvar $Di $Df\n";
         	   }
         	}
         }
         @Instance = ("WFG1", "WFG2", "WFG3", "WFG4", "WFG5", "WFG6", "WFG7", "WFG8", "WFG9");
            foreach(@Instance)
            {
            	for(my $nobj = 2; $nobj <=3; $nobj++)
            	{
            	   for(my $Sed = 1; $Sed <=35; $Sed++) ##Realizar 35 ejecuciones con distintas semilla de cada instancia..
            	   {
            	   	print $fout "~$PathAlgorithm/Ejecutable $PathAlgorithm $_ $Sed $nobj $pops $max_nfes $CR $F $nvar $Di $Df\n";
            	   }
            	}
            }
      @Instance = ("UF1", "UF2", "UF3", "UF4", "UF5", "UF6", "UF7");
      my $nobj=2;
         foreach(@Instance)
         {
         	   for(my $Sed = 1; $Sed <=35; $Sed++) ##Realizar 35 ejecuciones con distintas semilla de cada instancia..
         	   {
         	   	print $fout "~$PathAlgorithm/Ejecutable $PathAlgorithm $_ $Sed $nobj $pops $max_nfes $CR $F $nvar $Di $Df\n";
         	   }
         }
      
      @Instance = ("UF8", "UF9", "UF10");
      $nobj=3;
         foreach(@Instance)
         {
         	   for(my $Sed = 1; $Sed <=35; $Sed++) ##Realizar 35 ejecuciones con distintas semilla de cada instancia..
         	   {
         	   	print $fout "~$PathAlgorithm/Ejecutable $PathAlgorithm $_ $Sed $nobj $pops $max_nfes $CR $F $nvar $Di $Df\n";
         	   }
         }
}
