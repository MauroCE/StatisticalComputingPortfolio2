#!/bin/bash

for name in "setosa" "versicolor" "virginica"
do
  echo $name
  for col in {1..4}
  do
    cat iris.csv | grep -i $name | cut -d ',' -f $col | paste -sd+ - | bc | awk '{print $1/50}'
  done
done

