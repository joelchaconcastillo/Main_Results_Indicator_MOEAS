
PATH1=../POF
#PATH1=../POF_middle_eps
#PATH1=../POF_middle_Active_archive
DI=0.4
DF=0.5
for px in 0.2 # 0.4 0.6 0.8 1.0
do

  for instance in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF1 UF2 UF3 UF4 UF5 UF6 UF7 IMB1 IMB2 IMB3 IMB7 IMB8 IMB9
  do
     for sed in {1..35}
     do
     #tail -100 $PATH1/POF_*_${instance}_RUN${sed}_seed_*nobj_2_*_DI_${DI}*Px_${px}* | cut -f1,3 -d' '| awk '{ print sprintf("%.9f %.9f", $1, $2); }'   > POF/${instance}_2_${sed}_${px}
     #tail -100 $PATH1/POF_*_${instance}_RUN${sed}_seed_*nobj_2_*_DI_${DI}*DF_${DF}*Px_${px}* | cut -f1,3 -d' '  > POF/${instance}_2_${sed}_${px}
     tail -100 $PATH1/POF_*_${instance}_RUN${sed}_seed_*nobj_2_*_DI_${DI}*DF_${DF}*Px_${px}* | cut -f13,15 -d' '  > POF/${instance}_2_${sed}_${px}
     done
  done
  for instance in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF8 UF9 UF10 IMB4 IMB5 IMB6 IMB10
  do
     for sed in {1..35}
     do
     #tail -100 $PATH1/POF_*_${instance}_RUN${sed}_seed_*nobj_3*_DI_${DI}*DF_${DF}*Px_${px}* | cut -f1,3,5 -d' ' | awk '{ print sprintf("%.9f %.9f %.9f", $1, $2, $3); }'  > POF/${instance}_3_${sed}_${px}
     #tail -100 $PATH1/POF_*_${instance}_RUN${sed}_seed_*nobj_3*_DI_${DI}*DF_${DF}*Px_${px}* | cut -f1,3,5 -d' '   > POF/${instance}_3_${sed}_${px}
     tail -100 $PATH1/POF_*_${instance}_RUN${sed}_seed_*nobj_3*_DI_${DI}*DF_${DF}*Px_${px}* | cut -f19,21,23 -d' '  > POF/${instance}_3_${sed}_${px}
     done
  done
done
