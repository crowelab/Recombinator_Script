###############################
# HERE IS THE SCRIPT FOR IGHJ #
###############################
gawk  ' BEGIN{ RS = ">"; 
               ORS = "";
         }
 $1~/IGHJ/{ 
        if (match($1,/IGHJ[0-9]+\/[A-Z]+[0-9]*\-[0-9]+/)){ 
            String=substr($1,RSTART,RLENGTH)
            gsub(/[\/]/,"",String)
            #printf("%s\t%s\n",$1,String)
        }else if (match($1,/IGHJ[0-9]+\-[0-9]+\-[0-9]+/)){
            String=substr($1,RSTART,RLENGTH)
            #printf("%s\t%s\n",$1,String)
        }else if (match($1,/IGHJ[0-9]+\-[A-Z]*[0-9]+[A-Z]*/)){
            String=substr($1,RSTART,RLENGTH)
            #printf("%s\t%s\n",$1,String)
        }else if (match($1,/IGHJ[0-9]+/)){
            String=substr($1,RSTART,RLENGTH)
        }
        Sequence="";
        for(i=2;i<=NF;i++){
          Sequence=Sequence""$i
         }
        print $1" "length(Sequence)"\n"
        #printf("%s\n",String)
        Command=sprintf(" if [ ! -d %s ]; then mkdir %s; fi ",String,String);
        system(Command);
        Outfile=sprintf("%s/%s.fasta",String,String);
        gsub(/[\/]/,"",$1)
        printf(">%s\n%s\n",$1,Sequence) >> Outfile
      }' $1
