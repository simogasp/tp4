#!/usr/bin/env bash

declare -a arr=("checkout" "commit" "merge")

## now loop through the above array
for i in "${arr[@]}"
do
   echo "$i"
   filename="post-${i}"
   cp post-xxx-sample.txt ../.git/hooks/"${filename}"
   chmod +x ../.git/hooks/"${filename}"
done