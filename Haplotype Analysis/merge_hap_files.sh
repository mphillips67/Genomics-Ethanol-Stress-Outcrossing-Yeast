
#exaple of how to do this for High week 7 hap esitmates once they are generated. Can modify for any treatment cycle combination once all files are generated. 

head -n1 alt_ETH_rep49_wk07_C01_hap_freq.txt > High_wk7_haps.txt

for fname in *wk07*.txt
do
    tail -n+2 $fname >> High_wk7_haps.txt
done
