#!/bin/bash
again=yes 
while [[ $again == "yes" ]]
do
  echo "enter your name:"
  read name
  if [[ $name == "" ]]
  then
    break
  fi
  echo "enter your age:"
  read age
  if [[ $age -eq 0 ]]
  then
    break
  elif [[ $age -le 16 ]]
  then
    echo "$name,  your group is child"
  elif [[ $age -le 25 ]]
  then
    echo "$name,  your group is youth"
  еlse
    echo "$name,  your group is adult"
  fi
done
