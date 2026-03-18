#!/bin/bash

if [ -e clear_hanging.lock ]; then
    exit
fi

echo 1 > clear_hanging.lock

sleep 600
while true; do

    onemin=$(date --date="20 seconds ago" +%s);
    nrunning=$(ps aux | grep gulls_general | wc -l)
    date >> clear_hanging.lock
    if [ $nrunning -eq 1 ]; then
	rm clear_hanging.lock
	exit
    fi
    anykilled=0
    ls *.log | while read line; do
	ftime=$(date -r "$line" +%s);
	echo $onemin $ftime >> clear_hanging.lock;
	if [ $ftime -le $onemin ]; then
	    echo $line old >> clear_hanging.lock;
	    numtmp=${line%.log};
	    num=${numtmp##*_};
	    echo $num >> clear_hanging.lock;
	    pid=$(ps aux | grep gulls_general | grep " "$num'$' | awk '{print $2}');
	    echo kill $num $pid >> clear_hanging.lock;
	    kill -9 $pid;
	    anykilled=$((anykilled+1))
	fi
    done
    if [ $anykilled -gt 0 ]; then
	sleep 600
    else
	sleep 25
    fi
done
