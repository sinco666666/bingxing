# example.sh
#!/bin/sh
# PBS -N arm

pssh -h $PBS_NODEFILE mkdir -p /home/s2213628 1>&2
scp master:/home/s2213628/arm /home/s2213628
pscp -h $PBS_NODEFILE master:/home/s2213628/arm /home/s2213628 1>&2
/home/s2213628/arm
