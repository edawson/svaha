svaha - generate variation graphs for structural variants.
-----------
Make variation graphs from structural variants:  
[x] Deletions
[x] Inversions
[x] Insertions
[x] SNPs
[ ] Duplications
[ ] Transversions
[ ] Breakpoints


Don't worry: we'll be adding these as time permits.  


## What is svaha?
svaha is a small program that converts Variant Call Format (VCF) records into [Graphical Fragment Assembly](https://github.com/pmelsted/GFA-spec) format (i.e. sequence graphs like those in [vg](https://github.com/vgteam/vg)). It does so using a minimal single-base graph representation, the world's smallest and least-safe VCF parser (well, probably), and almost no dependencies.

## Build it
svaha brings in its own libraries, except for zlib. Make sure to have zlib installed.
It uses a frozen version of htslib and floating versions of gfakluge. To build svaha:  

                git clone --recursive https://github.com/edawson/svaha
                make


and that should do it.

## Run svaha
svaha takes a FASTA file and a VCF as arguments:  
        ```./svaha -r MYFASTA.fa -v MYVARIATION.vcf```

and outputs sorted GFA, which is text-based and easily exchangeable to other, more useful programs (like vg).


## Options
`-r`: a fasta reference  
`-v`: a vcf containing variants (must be relative to the given fasta)  
`-m`: maximum node size. When creating graphs for vg, make sure to use a maximum node size of between 32 and 1023.  
1023 is a hard limit (nothing 1024 or over will be indexable) and below 32 the graph begins to eat tons of memory. I tend to use `-m 64` or `-m 128`.

## Workflows

1. Build a variation graph with svaha containing structural variants  
2. Reduce node size with a ``cat result.gfa | vg view -F -v - | vg mod -X 1000 - > new_graph.vg`` to make the resulting graph indexable with GCSA2.  
3. Map reads to that graph using ```vg map``` 
4. Call variants using ```vg call``` or ```vg genotype```


## Get help
Reach out to me (@edawson) on GitHub and I'll do my best to help!
