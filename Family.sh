 gawk  ' BEGIN{ RS = ">"; 
               ORS = "";
         }
 $1~/IGHV/{ 
        if (match($1,/IGHV[1-9]+\//)){ 
            String=substr($1,RSTART,RLENGTH)
            gsub(/[\/]/,"",String)
        }else if (match($1,/IGHV[1-9]+/)){
            String=substr($1,RSTART,RLENGTH)
        }
        Sequence="";
        for(i=2;i<=NF;i++){
          Sequence=Sequence""$i
         }
      
        Command=sprintf(" if [ ! -d %s ]; then mkdir %s; fi ",String,String);
        system(Command);
        Outfile=sprintf("%s/%s.fasta",String,String);
        gsub(/[\/]/,"",$1)
        print $1" "length(Sequence)"\n"
        printf(">%s\n%s\n",$1,Sequence) >> Outfile
      }' $1

