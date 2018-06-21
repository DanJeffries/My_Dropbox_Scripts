#!bin/bash

name=$1
string1="Wake up, $name"
string2="The Matrix has you ..."
string3="Follow the white rabbit."
string4="Knock, knock, $name"

clear

for i in $string1
do
for j in $(echo $i | grep -o .)
do
echo -n  $j
sleep $(python -c "import random;print(random.uniform(0.2, 0.6))")
done
printf " "
done


sleep 3

clear

for i in $string2
do
for j in $(echo $i | grep -o .)
do
echo -n  $j
sleep $(python -c "import random;print(random.uniform(0.2, 0.6))")
done
printf " "
done


sleep 3
clear
         
for i in $string3
do
for j in $(echo $i | grep -o .)
do
echo -n  $j
sleep $(python -c "import random;print(random.uniform(0.2, 0.6))")
done
printf " "
done

sleep 4

clear

echo $string4
sleep 6         
