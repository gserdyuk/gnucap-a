#!/bin/bash

# calling <mydiff file dir> or <compare file file>
function mydiff {
    first=$1
    second=$2
    diff -u --ignore-matching-lines="Gnucap gnucap" $first $second 
    return $?
    }

#calling <printresult topic code>
function printresult {
   topic=$1
   code=$2
# o - ok, !=0 - trouble
   echo "$topic" "$code" 
   }
   
  