for i in `dir phispy/`; do file=`readlink -f phispy/$i/*_PhiSpy-prophage-coordinates.tsv`; count=`cat $file | grep -v '#prophage' | wc -l`; echo $i,$count; done  | sed 's/GCF/GCA/' > file

## In Saureus it would be necessary to do
cat file | sed 's/GCA_000953255.1/GCF_000953255.1/' | sed 's/GCA_000160335/GCF_000160335.2/' | file2
