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
5. Redesign the part relating to .hashmap files and how we iterate through a new input file
    - *Done*
    - It isn't faster, but it's waaay more memory-efficient.
    - A run on a single 8.6Gb (gzipped) file used to take 54Gb.
    - Now every run, regardless of input file size, uses at most around 5-6Gb RAM.
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