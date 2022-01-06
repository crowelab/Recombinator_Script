#gawk -F, '$1~/IGHV/ && $(NF-3)~/\*/{ 
gawk -F, '$1~/IGHV/{ 
           #print $1" "$(NF-3)
           #print $(NF-4)
           M=split($(NF-4),array,"")
           N=split($(NF-3),array2,"")
           Counter=1;
           for(i=1;i<=M;i=i+3){
             CODON[Counter]=array[i+0]""array[i+1]""array[i+2];
             Counter=Counter+1; 
           }
           for(j=1;j<Counter;j++){
             #print CODON[j]" "array2[j] 
           }
           Found=0;
           for(j=Counter;j>=1;j--){
              if (CODON[j]=="TGT"){
                  #print $1" "$(NF-3)
                  Found=1
                  j=1;continue
               }
           }
           if (Found==0){
               print $1 " TERMINAL TGT NOT FOUND"
           }
           delete array
           delete array2
           #print "----------------------------------------------------------------"
    }' Heavy.csv
