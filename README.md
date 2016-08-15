Lasso
-----------
Make variation graphs from structural variants:  
[x] Deletions
[ ] Inversions
[ ] Duplications
[ ] Transversions
[ ] Insertions
[ ] Breakpoints

Don't worry: we'll be adding these as time permits.  


## What is lasso?
lasso is a small program that converts Variant Call Format (VCF) records into [Graphical Fragment Assembly](https://github.com/pmelsted/GFA-spec) format (i.e. sequence graphs like those in [vg](https://github.com/vgteam/vg)). It does so using a minimal single-base graph representation, the world's smallest and least-safe VCF parser (well, probably), and almost no dependencies.

## Build it
lasso brings in its own libraries, except for zlib. Make sure to have zlib installed.
It uses a frozen version of htslib and floating versions of gfakluge. To build lasso:  

                git clone --recursive https://github.com/edawson/lasso
                make


and that should do it.

## Run lasso
lasso takes a FASTA file and a VCF as arguments:  
        ```./lasso -r MYFASTA.fa -v MYVARIATION.vcf```

and outputs sorted GFA, which is text-based and easily exchangeable to other, more useful programs (like vg).

## Workflows

1. Build a variation graph with lasso containing structural variants  
2. Reduce node size with a ``cat result.gfa | vg mod -X 1000 - > new_graph.vg`` to make the resulting graph indexable with GCSA2.  
3. Map reads to that graph using ```vg map``` 
4. Call variants using ```vg call``` or ```vg genotype```
