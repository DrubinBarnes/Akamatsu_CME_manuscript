#!/bin/bash
python ../py/preconfig "$1.cym.tpl"
mkdir ../simulations/$1
mv "$1"0* ../simulations/$1
echo "made and moved $1"
read -n1 -p "Do you want to send it up to the server? [Y/N] " answer
case $answer in
	Y | y) echo
		   echo "sending up to server"
		   
		   scp -r ../simulations/$1 mattakamatsu@dtn.brc.berkeley.edu:/global/scratch/mattakamatsu/;;

	N | n) echo 
		   echo "goodbye"
	exit;;
esac

read -n1 -p "Do you want to visit the server? [Y/N] " answer
case $answer in
	Y | y) echo
		   echo "visiting server"
		   
		   ssh mattakamatsu@hpc.brc.berkeley.edu;;

	N | n) echo 
		   echo "goodbye"
	exit;;
esac

# this part doesn't show up once youre on the server
# read -n1 -p "Do you want to execute the simulation? [Y/N]" answer
# case $answer in
# 	Y | y) echo
# 		   read -p "type how many repeats of each cym file (ideally number of cym files x repeats is a multiple of 24)" answer2
# 		   ./savioParallelCym.sh $1 $answer2;;
		   

# 	N | n) echo 
# 		   echo "goodbye"
# 	exit;;
# esac

