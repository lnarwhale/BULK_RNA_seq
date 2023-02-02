#!/bin/bash

function re_seq(){
	#--function:get the '-'sequence according to '+'sequence
	#--input:'+'sequence
	#--output:'-'sequence
	#--need:1.input sequence
	sequence=$1
	tmp1=$(echo ${sequence//A/t})
	tmp2=$(echo ${tmp1//C/g})
	tmp3=$(echo ${tmp2//G/c})
	tmp4=$(echo ${tmp3//T/a})
	tmp21=$(echo ${tmp4//a/A})
	tmp22=$(echo ${tmp21//c/C})
	tmp23=$(echo ${tmp22//t/T})
	tmp24=$(echo ${tmp23//g/G})
	echo $tmp24|rev
}
