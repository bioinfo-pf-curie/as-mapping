#!/bin/bash

snpsplitLog=$1

tot=$(grep 'total' ${snpsplitLog} | awk -F":" '{print $2}' | sed -e 's/\t//g')
unass=$(grep 'unassignable' ${snpsplitLog} | awk -F":" '{print $2}' | awk -F"(" '{print $1}' | sed -e 's/\t//g')
g1=$(grep 'specific for genome 1'  ${snpsplitLog} | awk -F":" '{print $2}' | awk -F"(" '{print $1}' | sed -e 's/\t//g')
g2=$(grep 'specific for genome 2'  ${snpsplitLog} | awk -F":" '{print $2}' | awk -F"(" '{print $1}' | sed -e 's/\t//g')
conf=$(grep 'conflicting'  ${snpsplitLog} | awk -F":" '{print $2}' | awk -F"(" '{print $1}' | sed -e 's/\t//g')

echo -e "Unassigned,$unass"
echo -e "Genome 1,$g1"
echo -e "Genome 2,$g2"
echo -e "Conflicting,$conf"
