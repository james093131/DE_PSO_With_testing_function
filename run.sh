   # ./MAIN 1 10000 10 100 100000 B 0.6 0.1 D >> B10_OG_DE_NEW.txt

   # ./MAIN 1 10000 10 100 100000 R 0.6 0.1 D >> R10_OG_DE_NEW.txt

   # ./MAIN 1 10000 10 100 100000 RO 0.6 0.1 D >> RO10_OG_DE_NEW.txt

   # ./MAIN 1 10000 10 100 100000 Z 0.6 0.1 D >>  Z10_OG_DE_NEW.txt


   j=30
   echo $i
   for (( j; j<31; j=j+1 ))
   do
      echo $j
      ./MAIN 10 4000 10 50 100000 $j 0.1 0.7 D >>  CEC_2017/10D/DE_CEC_${j}_T1.txt
      # ./MAIN 30 2000 10 100 100000 $j 0.7 0.2 P >>  CEC_test/PSO_CEC_${j}.txt
   done
 
