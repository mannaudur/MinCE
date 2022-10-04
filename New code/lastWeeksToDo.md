## Setting up MinCE 207

1. Set up MinCE on Tölva
    - *Done*
2. Start sketching all members of MinCE_207
    - k=31, s=5000 ?
3. Think about better ways to cluster->clique those members
    - Use time spent crunching on other things
    - Is it possible to use clustering algorithms on this?
    - Answer might be found by turning the problem on its head
4. Rewrite Sequences.cpp to run faster, it's much slower than it has to be.
    - Choosing fasta-sequencess can be moved outside of general scope
    - Simplify and classify
    - Sequences not working with batch file input
        - trying: for i in `cat helper_textfiles/cliques_all.txt`; do ../Code_MinCE/local/bin/sequences -N 20 ${i} ; done
5. Redesign the part relating to .hashmap files and how we iterate through a new input file
    - *Done*
    - It isn't faster, but it's waaay more memory-efficient.
    - A run on a single 8.6Gb (gzipped) file used to take 54Gb.
    - Now every run, regardless of input file size, uses at most around 5-6Gb RAM.
        - Except they don't, on Tölva it's taking 35 gigs, which is better than 54 before, by why?
        - I ran the same file three times on the MacBook and it only used 4.6, 5,5 and 6.1 Gb.
        - How come it varies so much? And why is it different on Tölva?
6. Write one final script, which takes in a results.tsv file and two parameters, T and R:
    - *Done*
7. Implement MinCE to work on many threads at once.
    - Don't know how to do it but it shouldn't be too difficult.

## Interpreting results from simulated data

1. Load kallisto index of genomes in results from mincing simulated data
    - For supposed 'true positives', both high- and low-ranking:
        - Check whether kallisto agrees that they're barely there
    - For high-ranking 'false positives':
        - Check whether I can align them appropriately on kallisto's index with the simulated data
    - Try to find some good threshold to use for future
2. Load simulated dataset into Kraken
    - I'll have to google that a bit more...
3. Try other simulated data?
    - Generate it myself based on wgsim?
    - Other datasets available?
4. Simulate other tests / calculate probability
    - How much mutation required for 5% chance of 4999, etc.


## Future

1. Setting up viral database 
    - already up on Tölva, hasn't been tested though
    - the one on Tölva is created from complete RefSeq alignments, there might be a better database out there
2. Expanding the database to include more bacterial genomes
    - Current one is based on GTDB v207 release
    - There is one purportedly containing over 660,000 bacterial genomes
        - Is it curated properly? How redundant is it?
        - https://www.sanger.ac.uk/news_item/new-database-of-660000-assembled-bacterial-genomes-sheds-light-on-the-evolution-of-bacteria/
        - https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001421
        - http://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k/

