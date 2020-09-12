#!/bin/bash

usage=$(cat << EOF

   # "This script calculates free binding energy for smallest and largest CAVER's tunnels"
   
   script.sh -p CAVER_PDB_DIR -i CAVER_OUT_DIR -o OUTNAME -c CPU -a C:N:O -t bottleneck -S sshlogin
          
      -p   : "CAVER's md_snapshots DIRECTORY"
      -i   : "CAVER's out DIRECTORY"
      -o   : "Output name"
      -c   : "Number of CPU (default: 1)"
      -a   : "Atoms types (default: C:N:O)"
      -t   : "curvature:length:bottleneck"
      -S   : "Distribute jobs to remote computers. The jobs will be run on a list of remote computers. (Example: name@nodeXX)"

    Atoms types:
	H      # Non H-bonding Hydrogen
	HD     # Donor 1 H-bond Hydrogen
	HS     # Donor S Spherical Hydrogen
	C      # Non H-bonding Aliphatic Carbon
	A      # Non H-bonding Aromatic Carbon
	N      # Non H-bonding Nitrogen
	NA     # Acceptor 1 H-bond Nitrogen
	NS     # Acceptor S Spherical Nitrogen
	OA     # Acceptor 2 H-bonds Oxygen
	OS     # Acceptor S Spherical Oxygen
	F      # Non H-bonding Fluorine
	Mg     # Non H-bonding Magnesium
	MG     # Non H-bonding Magnesium
	P      # Non H-bonding Phosphorus
	SA     # Acceptor 2 H-bonds Sulphur
	S      # Non H-bonding Sulphur
	Cl     # Non H-bonding Chlorine
	CL     # Non H-bonding Chlorine
	Ca     # Non H-bonding Calcium
	CA     # Non H-bonding Calcium
	Mn     # Non H-bonding Manganese
	MN     # Non H-bonding Manganese
	Fe     # Non H-bonding Iron
	FE     # Non H-bonding Iron
	Zn     # Non H-bonding Zinc
	ZN     # Non H-bonding Zinc
	Br     # Non H-bonding Bromine
	BR     # Non H-bonding Bromine
	I      # Non H-bonding Iodine
	Z      # Non H-bonding covalent map
	G      # Ring closure Glue Aliphatic Carbon  # SF
	GA     # Ring closure Glue Aromatic Carbon   # SF
	J      # Ring closure Glue Aliphatic Carbon  # SF
	Q      # Ring closure Glue Aliphatic Carbon  # SF

    Option:
 
      -h : help message

    Version:

      version : 1.0
      author : Javier Rodriguez-Salarichs

EOF
);


while getopts ":p:i:o:c:a:t:S:h" opt; do
  case $opt in
    p)
      INPUT=$OPTARG
      ;;
    i)
      IOUT=$OPTARG
      ;;
    o)
      OUTPUT=$OPTARG
      ;;
    c)
      CPU=$OPTARG
      ;;
    a)
      AT_ALL=$OPTARG
      ;;
    t)
      bottleneck_or_length=$OPTARG
      ;;
    S)
      remote_node=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      echo "$usage"
      exit 1;
      ;;
    h)
       echo "$usage"
       exit 1;
       ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1;
      ;;
  esac
done

if [ -z $AT_ALL ];then
AT_ALL="C:N:O"
fi

if [ -z $CPU ];then
CPU=1
fi

if [ -z $INPUT ];then
echo -e "\tIt needs to define the CAVER's md_snapshots DIRECTORY.\n\tType -h for more options." >&2
      exit 1;
fi

if [ -z $IOUT ];then
echo -e "\tIt needs to define the CAVER's out DIRECTORY.\n\tType -h for more options." >&2
      exit 1;
fi

if [ ! -d "$INPUT" ];then
  echo -e "\tout DIRECTORY doesn't exist.\n\tType -h for more options.">&2
  exit 1
fi

if [ ! -d "$IOUT" ];then
  echo -e "\tout DIRECTORY doesn't exist.\n\tType -h for more options.">&2
  exit 1
fi

if [ -z $OUTPUT ];then
echo -e "\tIt needs to define the output file.\n\tType -h for more options." >&2
      exit 1;
fi

if [ -a $OUTPUT ];then
echo "Output file exists." >&2
      exit 1;
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

T="$(date +%s)"

main_funct () {

dir_nobin=`echo $DIR|awk '{sub("/bin","",$0);print $0}'`

mkdir $OUTPUT
mkdir $OUTPUT"/out"

mkdir $OUTPUT"/inputs"
mkdir $OUTPUT"/inputs/protein"
mkdir $OUTPUT"/inputs/tunnels"

T=`echo $AT_ALL |awk -F":" '{print NF}'`
i=1
while [ $i -le $T ]
do
 AT[$i]=`echo $AT_ALL |awk -v n=$i -F":" '{print $n }'`
 mkdir $OUTPUT"/out/"${AT[$i]}"_file"
i=$[$i+1]
done

mkdir $OUTPUT"/tmp"

#echo $DIR
#for i in `ls $INPUT/`;do
# echo "$i: from pdb to pdbqt"
# $DIR"/prepare_receptor-ipdb-opdbqt.py" -r $INPUT"/"$i -o $OUTPUT"/inputs/"$i"qt"
#done

ls $INPUT/ &>$OUTPUT/tmp/int$$.txt 

$dir_nobin"/third-party/parallel-20170622/bin/parallel" -k -j"$CPU" $DIR"/prepare_receptor-ipdb-opdbqt.py" -r $INPUT"/"{} -o $OUTPUT"/inputs/protein/"{}"qt" < $OUTPUT/tmp/int$$.txt

rm $OUTPUT/tmp/int$$.txt

Tl=`ls $IOUT/data/clusters/ |wc -l`

l=1
for j in `ls $IOUT/data/clusters/|awk '{sub(".pdb","",$0);print $0}'`
do
 mkdir $OUTPUT"/inputs/tunnels/tunnel_"$l
 i=1
 while [ $i -le $T ]
 do

  mkdir $OUTPUT"/inputs/tunnels/tunnel_"$l"/"${AT[$i]}"_file"

#  cat $IOUT"/data/clusters/"$j".pdb" | parallel --pipe -k -j"$CPU" awk -v AT=${AT[$i]} -v NAME=$j -v OUTP=$OUTPUT"/inputs/tunnels/tunnel_"$l"/"${AT[$i]}"_file" \''{if(/ATOM/){sub("H",AT,$3);N+=1;printf "%s\n%s %6s %2s %5s %4s %12.3f %7.3f %7.3f %5s %5s %9s %s\n%s\n%s\n","ROOT","ATOM","1",AT,"UNK","1",$7,$8,$9,"0.00","0.00","+0.000",AT,"ENDROOT","TORSDOF 0" > OUTP"/"NAME"_"AT"_"N".pdbqt" };if(length==0){exit}}'\'  

 i=$[$i+1]
 done
l=$[$l+1]
done


col_rl=6
if [ -z $bottleneck_or_length ];then
 bottleneck_or_length="bottleneck"
fi

bottleneck_or_length=`echo $bottleneck_or_length | awk '{ print tolower($0) }'`

if [ $bottleneck_or_length == "length" ];then
 col_rl=8
elif [ $bottleneck_or_length == "bottleneck" ];then 
 col_rl=6
elif [ $bottleneck_or_length == "Curvature" ];then  
 col_rl=9
fi

i=1
while [ $i -le $T ]
do
 awk -v col_rl=$col_rl -v AT=${AT[$i]} -v OUTP=$OUTPUT"/inputs/tunnels" -v OUTT=$OUTPUT"/tmp" 'BEGIN{FS=",";cluster=0}{if(NR==FNR){if(cluster!=$2){cluster=int($2);min[cluster,2]=100000};if(max[cluster,2]<$col_rl){max[cluster,2]=$col_rl;max[cluster,1]=$1};if(min[cluster,2]>$col_rl){min[cluster,2]=$col_rl;min[cluster,1]=$1};next};if(NR!=FNR){for(i=1;i<=cluster;i++){if(min[i,1]==$1 && i==$2){   if($13==" R"){Ts[i]=NF;for(j=14;j<=NF;j++){Rs[i,j]=$j}};if($13==" X"){for(j=14;j<=NF;j++){Xs[i,j]=$j}};if($13==" Y"){for(j=14;j<=NF;j++){Ys[i,j]=$j}};if($13==" Z"){for(j=14;j<=NF;j++){Zs[i,j]=$j}};if($13==" distance"){for(j=14;j<=NF;j++){dists[i,j]=$j}};if($13==" length"){for(j=14;j<=NF;j++){lens[i,j]=$j}}    };if(max[i,1]==$1 && i==$2){   if($13==" R"){Tl[i]=NF;for(j=14;j<=NF;j++){Rl[i,j]=$j}};if($13==" X"){for(j=14;j<=NF;j++){Xl[i,j]=$j}};if($13==" Y"){for(j=14;j<=NF;j++){Yl[i,j]=$j}};if($13==" Z"){for(j=14;j<=NF;j++){Zl[i,j]=$j}};if($13==" distance"){for(j=14;j<=NF;j++){distl[i,j]=$j}};if($13==" length"){for(j=14;j<=NF;j++){lenl[i,j]=$j}}    }}};next}END{ for(i=1;i<=cluster;i++){print i,min[i,1],min[i,2]"\n Point Length Radii Curvature" > OUTT"/tunnel_smallest_"i".txt" ;print i,max[i,1],max[i,2]"\n Point Length Radii Curvature" > OUTT"/tunnel_largest_"i".txt" ;for(j=14;j<=Ts[i];j++){ printf "%s\n%s %6s %2s %5s %4s %12.3f %7.3f %7.3f %5s %5s %9s %s\n%s\n%s\n","ROOT","ATOM","1",AT,"UNK","1",Xs[i,j],Ys[i,j],Zs[i,j],"0.00","0.00","+0.000",AT,"ENDROOT","TORSDOF 0" > OUTP"/tunnel_"i"/"AT"_file/tunnel_smallest_"i"_"AT"_"j-13".pdbqt" ; print j-13,lens[i,j],Rs[i,j],lens[i,j]/dists[i,j] > OUTT"/tunnel_smallest_"i".txt" };for(j=14;j<=Tl[i];j++){printf "%s\n%s %6s %2s %5s %4s %12.3f %7.3f %7.3f %5s %5s %9s %s\n%s\n%s\n","ROOT","ATOM","1",AT,"UNK","1",Xl[i,j],Yl[i,j],Zl[i,j],"0.00","0.00","+0.000",AT,"ENDROOT","TORSDOF 0" > OUTP"/tunnel_"i"/"AT"_file/tunnel_largest_"i"_"AT"_"j-13".pdbqt";print j-13,lenl[i,j],Rl[i,j],lenl[i,j]/distl[i,j] > OUTT"/tunnel_largest_"i".txt"  } }}' $IOUT"/analysis/tunnel_characteristics.csv"  $IOUT"/analysis/tunnel_profiles.csv"

 i=$[$i+1]
done


Tt=`echo $l|awk '{print $l-1}'`

mkdir $OUTPUT"/out/pdb"

for i in `ls $OUTPUT"/inputs/tunnels/"`
do
 mkdir $OUTPUT"/tmp/"$i
# mkdir $OUTPUT"/out/"$i

 kk=`echo $i | awk '{sub("tunnel_","",$0);print $0}'`

 for j in `ls $OUTPUT"/inputs/tunnels/"$i"/"`
 do
 
 mkdir $OUTPUT"/tmp/"$i"/"$j
# mkdir $OUTPUT"/out/"$i"/"$j

 $dir_nobin"/third-party/parallel-20170622/bin/parallel" -k -j"$CPU" -a <(ls $OUTPUT"/inputs/tunnels/"$i"/"$j"/")   $DIR"/vina" --cpu $CPU --score_only --receptor $OUTPUT"/inputs/protein/"{2} --ligand $OUTPUT"/inputs/tunnels/"$i"/"$j"/"{1}\|awk -v val1={1} -v val2={2} -v  OUTF=$OUTPUT"/tmp/"$i"/"$j"/"{1.}".txt" \''/Affinity/{print val1,val2,$2 >> OUTF }'\' :::: <(ls $OUTPUT"/inputs/protein/")

 $dir_nobin"/third-party/parallel-20170622/bin/parallel" -k -j"$CPU" -a <( ls $OUTPUT/tmp/$i/$j/*largest* ) awk -v OUTF=$OUTPUT"/tmp/"$i"/"$j -v tunnel=$i -v AT=$j -v name={1/.} \''{E[NR]=$3;sum+=$3}END{mean=sum/NR;sub("_file","",AT);for(i=1;i<=NR;i++){ sd+=((E[i]-mean)^2)};printf "%s %0.3f %s %0.3f\n",name,mean,"+/-",sqrt(sd/NR) >> OUTF"/tunnel_largest_"AT".txt" }'\' {1}

$dir_nobin"/third-party/parallel-20170622/bin/parallel" -k -j"$CPU" -a <( ls $OUTPUT/tmp/$i/$j/*smallest* ) awk -v OUTF=$OUTPUT"/tmp/"$i"/"$j -v tunnel=$i -v AT=$j -v name={1/.} \''{E[NR]=$3;sum+=$3}END{mean=sum/NR;sub("_file","",AT);for(i=1;i<=NR;i++){ sd+=((E[i]-mean)^2)};printf "%s %0.3f %s %0.3f\n",name,mean,"+/-",sqrt(sd/NR) >> OUTF"/tunnel_smallest_"AT".txt" }'\' {1}

 rr=`echo $j | awk '{sub("_file","",$0);print $0}'`

$dir_nobin"/third-party/parallel-20170622/bin/parallel" -k -j"$CPU" -a <( ls $OUTPUT/inputs/tunnels/$i/$j/*smallest* ) awk -v OUTF=$OUTPUT"/tmp/"$i"/"$j -v tunnel=$i -v AT=$rr -v name={1/.} \''FNR==NR{k+=1;name1[k]=$1;for(i=1;i<=length($1);i++){if(substr($1,i,1) == "_"){n=i}};num[k]=substr($1,n+1,length($1));E[k]=$2;Err[k]=$4;next}FNR!=NR{for(i=1;i<=k;i++){if(name == name1[i] && /ATOM/){ print num[i],tunnel,AT,name,$6,$7,$8,E[i],Err[i],exp(E[i]/0.593)}};next}'\' {2} {1} :::: <(ls $OUTPUT"/tmp/"$i"/"$j"/tunnel_smallest_"$rr".txt") |sort -g -k 1 > $OUTPUT"/tmp/"$i"/"$j"/"$i"_smallest_"$rr"_.txt"

$dir_nobin"/third-party/parallel-20170622/bin/parallel" -k -j"$CPU" -a <( ls $OUTPUT/inputs/tunnels/$i/$j/*largest* ) awk -v OUTF=$OUTPUT"/tmp/"$i"/"$j -v tunnel=$i -v AT=$rr -v name={1/.} \''FNR==NR{k+=1;name1[k]=$1;for(i=1;i<=length($1);i++){if(substr($1,i,1) == "_"){n=i}};num[k]=substr($1,n+1,length($1));E[k]=$2;Err[k]=$4;next}FNR!=NR{for(i=1;i<=k;i++){if(name == name1[i] && /ATOM/){ print num[i],tunnel,AT,name,$6,$7,$8,E[i],Err[i],exp(E[i]/0.593)}};next}'\' {2} {1} :::: <(ls $OUTPUT"/tmp/"$i"/"$j"/tunnel_largest_"$rr".txt") |sort -g -k 1 > $OUTPUT"/tmp/"$i"/"$j"/"$i"_largest_"$rr"_.txt"



awk '{printf "%-s%+5s  %-s %2s%3s%2s%4s%s   %8.3f%8.3f%8.3f%6.2f%6.2f%12s\n","ATOM  ",NR,"H"," ","FIL","T","1"," ",$5,$6,$7,"1.00",$8,"H"}END{for(i=2;i<=NR;i++){printf "%s%+5s%+5s\n","CONECT",i-1,i};print "END"}' $OUTPUT"/tmp/"$i"/"$j"/"$i"_smallest_"$rr"_.txt" > $OUTPUT"/out/pdb/"$i"_smallest_"$rr"_.pdb"

awk '{printf "%-s%+5s  %-s %2s%3s%2s%4s%s   %8.3f%8.3f%8.3f%6.2f%6.2f%12s\n","ATOM  ",NR,"H"," ","FIL","T","1"," ",$5,$6,$7,"1.00",$8,"H"}END{for(i=2;i<=NR;i++){printf "%s%+5s%+5s\n","CONECT",i-1,i};print "END"}' $OUTPUT"/tmp/"$i"/"$j"/"$i"_largest_"$rr"_.txt" > $OUTPUT"/out/pdb/"$i"_largest_"$rr"_.pdb"


awk -v OUTF=$OUTPUT"/out/"$j"/"$i"_smallest.txt" 'BEGIN{print "length\t radii[Angstrom]\t curvature\t DGb[kcal/mol]\t errDGb\t Kd" > OUTF }{if(NR==FNR){if(FNR>2){var[FNR-2,1]=$2;var[FNR-2,2]=$3;var[FNR-2,3]=$4};next};if(NR!=FNR){var[FNR,4]=$8;var[FNR,5]=$9;var[FNR,6]=$10;next}} END {for(i=1;i<=FNR;i++){print var[i,1]"\t"var[i,2]"\t"var[i,3]"\t"var[i,4]"\t"var[i,5]"\t"var[i,6] >> OUTF }}' $OUTPUT"/tmp/tunnel_smallest_"$kk".txt" $OUTPUT"/tmp/"$i"/"$j"/"$i"_smallest_"$rr"_.txt"

awk -v OUTF=$OUTPUT"/out/"$j"/"$i"_largest.txt" 'BEGIN{print "length\t radii[Angstrom]\t curvature\t DGb[kcal/mol]\t errDGb\t Kd" > OUTF }{if(NR==FNR){if(FNR>2){var[FNR-2,1]=$2;var[FNR-2,2]=$3;var[FNR-2,3]=$4};next};if(NR!=FNR){var[FNR,4]=$8;var[FNR,5]=$9;var[FNR,6]=$10;next}} END {for(i=1;i<=FNR;i++){print var[i,1]"\t"var[i,2]"\t"var[i,3]"\t"var[i,4]"\t"var[i,5]"\t"var[i,6] >> OUTF }}' $OUTPUT"/tmp/tunnel_largest_"$kk".txt" $OUTPUT"/tmp/"$i"/"$j"/"$i"_largest_"$rr"_.txt"


 done
done

}

if [ -z $remote_node ];then
main_funct
else
ssh -T -oLogLevel=quiet $remote_node "bash -s $INPUT $IOUT $OUTPUT $CPU $AT_ALL $bottleneck_or_length $DIR" << 'EOF'
INPUT=$1
IOUT=$2
OUTPUT=$3
CPU=$4
AT_ALL=$5
bottleneck_or_length=$6

DIR=$7
dir_nobin=`echo $DIR|awk '{sub("/bin","",$0);print $0}'`

mkdir $OUTPUT
mkdir $OUTPUT"/out"

mkdir $OUTPUT"/inputs"
mkdir $OUTPUT"/inputs/protein"
mkdir $OUTPUT"/inputs/tunnels"

T=`echo $AT_ALL |awk -F":" '{print NF}'`
i=1
while [ $i -le $T ]
do
 AT[$i]=`echo $AT_ALL |awk -v n=$i -F":" '{print $n }'`
 mkdir $OUTPUT"/out/"${AT[$i]}"_file"
i=$[$i+1]
done

mkdir $OUTPUT"/tmp"

#echo $DIR
#for i in `ls $INPUT/`;do
# echo "$i: from pdb to pdbqt"
# $DIR"/prepare_receptor-ipdb-opdbqt.py" -r $INPUT"/"$i -o $OUTPUT"/inputs/"$i"qt"
#done

ls $INPUT/ &>$OUTPUT/tmp/int$$.txt 

$dir_nobin"/third-party/parallel-20170622/bin/parallel" -k -j"$CPU" $DIR"/prepare_receptor-ipdb-opdbqt.py" -r $INPUT"/"{} -o $OUTPUT"/inputs/protein/"{}"qt" < $OUTPUT/tmp/int$$.txt

rm $OUTPUT/tmp/int$$.txt

Tl=`ls $IOUT/data/clusters/ |wc -l`

l=1
for j in `ls $IOUT/data/clusters/|awk '{sub(".pdb","",$0);print $0}'`
do
 mkdir $OUTPUT"/inputs/tunnels/tunnel_"$l
 i=1
 while [ $i -le $T ]
 do

  mkdir $OUTPUT"/inputs/tunnels/tunnel_"$l"/"${AT[$i]}"_file"

#  cat $IOUT"/data/clusters/"$j".pdb" | parallel --pipe -k -j"$CPU" awk -v AT=${AT[$i]} -v NAME=$j -v OUTP=$OUTPUT"/inputs/tunnels/tunnel_"$l"/"${AT[$i]}"_file" \''{if(/ATOM/){sub("H",AT,$3);N+=1;printf "%s\n%s %6s %2s %5s %4s %12.3f %7.3f %7.3f %5s %5s %9s %s\n%s\n%s\n","ROOT","ATOM","1",AT,"UNK","1",$7,$8,$9,"0.00","0.00","+0.000",AT,"ENDROOT","TORSDOF 0" > OUTP"/"NAME"_"AT"_"N".pdbqt" };if(length==0){exit}}'\'  

 i=$[$i+1]
 done
l=$[$l+1]
done


col_rl=6
if [ -z $bottleneck_or_length ];then
 bottleneck_or_length="bottleneck"
fi

bottleneck_or_length=`echo $bottleneck_or_length | awk '{ print tolower($0) }'`

if [ $bottleneck_or_length == "length" ];then
 col_rl=8
elif [ $bottleneck_or_length == "bottleneck" ];then 
 col_rl=6
elif [ $bottleneck_or_length == "Curvature" ];then  
 col_rl=9
fi

i=1
while [ $i -le $T ]
do
 awk -v col_rl=$col_rl -v AT=${AT[$i]} -v OUTP=$OUTPUT"/inputs/tunnels" -v OUTT=$OUTPUT"/tmp" 'BEGIN{FS=",";cluster=0}{if(NR==FNR){if(cluster!=$2){cluster=int($2);min[cluster,2]=100000};if(max[cluster,2]<$col_rl){max[cluster,2]=$col_rl;max[cluster,1]=$1};if(min[cluster,2]>$col_rl){min[cluster,2]=$col_rl;min[cluster,1]=$1};next};if(NR!=FNR){for(i=1;i<=cluster;i++){if(min[i,1]==$1 && i==$2){   if($13==" R"){Ts[i]=NF;for(j=14;j<=NF;j++){Rs[i,j]=$j}};if($13==" X"){for(j=14;j<=NF;j++){Xs[i,j]=$j}};if($13==" Y"){for(j=14;j<=NF;j++){Ys[i,j]=$j}};if($13==" Z"){for(j=14;j<=NF;j++){Zs[i,j]=$j}};if($13==" distance"){for(j=14;j<=NF;j++){dists[i,j]=$j}};if($13==" length"){for(j=14;j<=NF;j++){lens[i,j]=$j}}    };if(max[i,1]==$1 && i==$2){   if($13==" R"){Tl[i]=NF;for(j=14;j<=NF;j++){Rl[i,j]=$j}};if($13==" X"){for(j=14;j<=NF;j++){Xl[i,j]=$j}};if($13==" Y"){for(j=14;j<=NF;j++){Yl[i,j]=$j}};if($13==" Z"){for(j=14;j<=NF;j++){Zl[i,j]=$j}};if($13==" distance"){for(j=14;j<=NF;j++){distl[i,j]=$j}};if($13==" length"){for(j=14;j<=NF;j++){lenl[i,j]=$j}}    }}};next}END{ for(i=1;i<=cluster;i++){print i,min[i,1],min[i,2]"\n Point Length Radii Curvature" > OUTT"/tunnel_smallest_"i".txt" ;print i,max[i,1],max[i,2]"\n Point Length Radii Curvature" > OUTT"/tunnel_largest_"i".txt" ;for(j=14;j<=Ts[i];j++){ printf "%s\n%s %6s %2s %5s %4s %12.3f %7.3f %7.3f %5s %5s %9s %s\n%s\n%s\n","ROOT","ATOM","1",AT,"UNK","1",Xs[i,j],Ys[i,j],Zs[i,j],"0.00","0.00","+0.000",AT,"ENDROOT","TORSDOF 0" > OUTP"/tunnel_"i"/"AT"_file/tunnel_smallest_"i"_"AT"_"j-13".pdbqt" ; print j-13,lens[i,j],Rs[i,j],lens[i,j]/dists[i,j] > OUTT"/tunnel_smallest_"i".txt" };for(j=14;j<=Tl[i];j++){printf "%s\n%s %6s %2s %5s %4s %12.3f %7.3f %7.3f %5s %5s %9s %s\n%s\n%s\n","ROOT","ATOM","1",AT,"UNK","1",Xl[i,j],Yl[i,j],Zl[i,j],"0.00","0.00","+0.000",AT,"ENDROOT","TORSDOF 0" > OUTP"/tunnel_"i"/"AT"_file/tunnel_largest_"i"_"AT"_"j-13".pdbqt";print j-13,lenl[i,j],Rl[i,j],lenl[i,j]/distl[i,j] > OUTT"/tunnel_largest_"i".txt"  } }}' $IOUT"/analysis/tunnel_characteristics.csv"  $IOUT"/analysis/tunnel_profiles.csv"

 i=$[$i+1]
done


Tt=`echo $l|awk '{print $l-1}'`

mkdir $OUTPUT"/out/pdb"

for i in `ls $OUTPUT"/inputs/tunnels/"`
do
 mkdir $OUTPUT"/tmp/"$i
# mkdir $OUTPUT"/out/"$i

 kk=`echo $i | awk '{sub("tunnel_","",$0);print $0}'`

 for j in `ls $OUTPUT"/inputs/tunnels/"$i"/"`
 do
 
 mkdir $OUTPUT"/tmp/"$i"/"$j
# mkdir $OUTPUT"/out/"$i"/"$j

 $dir_nobin"/third-party/parallel-20170622/bin/parallel" -k -j"$CPU" -a <(ls $OUTPUT"/inputs/tunnels/"$i"/"$j"/")   $DIR"/vina" --cpu $CPU --score_only --receptor $OUTPUT"/inputs/protein/"{2} --ligand $OUTPUT"/inputs/tunnels/"$i"/"$j"/"{1}\|awk -v val1={1} -v val2={2} -v  OUTF=$OUTPUT"/tmp/"$i"/"$j"/"{1.}".txt" \''/Affinity/{print val1,val2,$2 >> OUTF }'\' :::: <(ls $OUTPUT"/inputs/protein/")

 $dir_nobin"/third-party/parallel-20170622/bin/parallel" -k -j"$CPU" -a <( ls $OUTPUT/tmp/$i/$j/*largest* ) awk -v OUTF=$OUTPUT"/tmp/"$i"/"$j -v tunnel=$i -v AT=$j -v name={1/.} \''{E[NR]=$3;sum+=$3}END{mean=sum/NR;sub("_file","",AT);for(i=1;i<=NR;i++){ sd+=((E[i]-mean)^2)};printf "%s %0.3f %s %0.3f\n",name,mean,"+/-",sqrt(sd/NR) >> OUTF"/tunnel_largest_"AT".txt" }'\' {1}

$dir_nobin"/third-party/parallel-20170622/bin/parallel" -k -j"$CPU" -a <( ls $OUTPUT/tmp/$i/$j/*smallest* ) awk -v OUTF=$OUTPUT"/tmp/"$i"/"$j -v tunnel=$i -v AT=$j -v name={1/.} \''{E[NR]=$3;sum+=$3}END{mean=sum/NR;sub("_file","",AT);for(i=1;i<=NR;i++){ sd+=((E[i]-mean)^2)};printf "%s %0.3f %s %0.3f\n",name,mean,"+/-",sqrt(sd/NR) >> OUTF"/tunnel_smallest_"AT".txt" }'\' {1}

 rr=`echo $j | awk '{sub("_file","",$0);print $0}'`

$dir_nobin"/third-party/parallel-20170622/bin/parallel" -k -j"$CPU" -a <( ls $OUTPUT/inputs/tunnels/$i/$j/*smallest* ) awk -v OUTF=$OUTPUT"/tmp/"$i"/"$j -v tunnel=$i -v AT=$rr -v name={1/.} \''FNR==NR{k+=1;name1[k]=$1;for(i=1;i<=length($1);i++){if(substr($1,i,1) == "_"){n=i}};num[k]=substr($1,n+1,length($1));E[k]=$2;Err[k]=$4;next}FNR!=NR{for(i=1;i<=k;i++){if(name == name1[i] && /ATOM/){ print num[i],tunnel,AT,name,$6,$7,$8,E[i],Err[i],exp(E[i]/0.593)}};next}'\' {2} {1} :::: <(ls $OUTPUT"/tmp/"$i"/"$j"/tunnel_smallest_"$rr".txt") |sort -g -k 1 > $OUTPUT"/tmp/"$i"/"$j"/"$i"_smallest_"$rr"_.txt"

$dir_nobin"/third-party/parallel-20170622/bin/parallel" -k -j"$CPU" -a <( ls $OUTPUT/inputs/tunnels/$i/$j/*largest* ) awk -v OUTF=$OUTPUT"/tmp/"$i"/"$j -v tunnel=$i -v AT=$rr -v name={1/.} \''FNR==NR{k+=1;name1[k]=$1;for(i=1;i<=length($1);i++){if(substr($1,i,1) == "_"){n=i}};num[k]=substr($1,n+1,length($1));E[k]=$2;Err[k]=$4;next}FNR!=NR{for(i=1;i<=k;i++){if(name == name1[i] && /ATOM/){ print num[i],tunnel,AT,name,$6,$7,$8,E[i],Err[i],exp(E[i]/0.593)}};next}'\' {2} {1} :::: <(ls $OUTPUT"/tmp/"$i"/"$j"/tunnel_largest_"$rr".txt") |sort -g -k 1 > $OUTPUT"/tmp/"$i"/"$j"/"$i"_largest_"$rr"_.txt"



awk '{printf "%-s%+5s  %-s %2s%3s%2s%4s%s   %8.3f%8.3f%8.3f%6.2f%6.2f%12s\n","ATOM  ",NR,"H"," ","FIL","T","1"," ",$5,$6,$7,"1.00",$8,"H"}END{for(i=2;i<=NR;i++){printf "%s%+5s%+5s\n","CONECT",i-1,i};print "END"}' $OUTPUT"/tmp/"$i"/"$j"/"$i"_smallest_"$rr"_.txt" > $OUTPUT"/out/pdb/"$i"_smallest_"$rr"_.pdb"

awk '{printf "%-s%+5s  %-s %2s%3s%2s%4s%s   %8.3f%8.3f%8.3f%6.2f%6.2f%12s\n","ATOM  ",NR,"H"," ","FIL","T","1"," ",$5,$6,$7,"1.00",$8,"H"}END{for(i=2;i<=NR;i++){printf "%s%+5s%+5s\n","CONECT",i-1,i};print "END"}' $OUTPUT"/tmp/"$i"/"$j"/"$i"_largest_"$rr"_.txt" > $OUTPUT"/out/pdb/"$i"_largest_"$rr"_.pdb"


awk -v OUTF=$OUTPUT"/out/"$j"/"$i"_smallest.txt" 'BEGIN{print "length\t radii[Angstrom]\t curvature\t DGb[kcal/mol]\t errDGb\t Kd" > OUTF }{if(NR==FNR){if(FNR>2){var[FNR-2,1]=$2;var[FNR-2,2]=$3;var[FNR-2,3]=$4};next};if(NR!=FNR){var[FNR,4]=$8;var[FNR,5]=$9;var[FNR,6]=$10;next}} END {for(i=1;i<=FNR;i++){print var[i,1]"\t"var[i,2]"\t"var[i,3]"\t"var[i,4]"\t"var[i,5]"\t"var[i,6] >> OUTF }}' $OUTPUT"/tmp/tunnel_smallest_"$kk".txt" $OUTPUT"/tmp/"$i"/"$j"/"$i"_smallest_"$rr"_.txt"

awk -v OUTF=$OUTPUT"/out/"$j"/"$i"_largest.txt" 'BEGIN{print "length\t radii[Angstrom]\t curvature\t DGb[kcal/mol]\t errDGb\t Kd" > OUTF }{if(NR==FNR){if(FNR>2){var[FNR-2,1]=$2;var[FNR-2,2]=$3;var[FNR-2,3]=$4};next};if(NR!=FNR){var[FNR,4]=$8;var[FNR,5]=$9;var[FNR,6]=$10;next}} END {for(i=1;i<=FNR;i++){print var[i,1]"\t"var[i,2]"\t"var[i,3]"\t"var[i,4]"\t"var[i,5]"\t"var[i,6] >> OUTF }}' $OUTPUT"/tmp/tunnel_largest_"$kk".txt" $OUTPUT"/tmp/"$i"/"$j"/"$i"_largest_"$rr"_.txt"


 done
done




EOF

fi

T="$(($(date +%s)-T))"

echo -e "\nTime in seconds: ${T}"

printf "Pretty format: %02d:%02d:%02d:%02d\n" "$((T/86400))" "$((T/3600%24))" "$((T/60%60))" "$((T%60))"

echo ""


exit

