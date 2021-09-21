PATHA=../HV/  #/home/joel.chacon/Current/MyResearchTopics/Data_Decomposition/StateOfTheArt/HV

   for i in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF1 UF2 UF3 UF4 UF5 UF6 UF7 IMB1 IMB2 IMB3 IMB7 IMB8 IMB9
   do
     echo "--"${i}_2; 
     for Di in 0.2 # 0.4 0.6 0.8 1.0
     do
       echo ${PATHA}/POF/${i}_2_${Di}
#       echo ${PATHA}/POF_sms_emoa/${i}_2_${Di}
     done
   done
   for i in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF8 UF9 UF10 IMB4 IMB5 IMB6 IMB10
   do
     echo "--"${i}_3; 
     for Di in 0.2 #0.4 0.6 0.8 1.0
     do
       echo ${PATHA}/POF/${i}_3_${Di}
 #      echo ${PATHA}/POF_sms_emoa/${i}_3_${Di}
     done
   done


