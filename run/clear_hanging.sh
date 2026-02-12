#!/bin/bash

fivemin=$(date --date="5 minutes ago" +%s);
ls *.log | while read line; do
    ftime=$(date -r "$line" +%s);
    echo $fivemin $ftime;
    if [ $ftime -le $fivemin ]; then
	echo $line old;
	numtmp=${line%.log};
	num=${numtmp##*_};
	echo $num;
	pid=$(ps aux | grep gulls_general | grep " "$num'$' | awk '{print $2}');
	echo $num $pid;
	kill -9 $pid;
    fi;
done
