## Setting up MinCE 207

1. Set up MinCE on Tölva
    - Once download is finished, do the same thing to Archaea
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
    - It has to be faster than this, though I don't immediately see any huge time-saves
6. Write one final script, which takes in a results.tsv file and two parameters, T and R:
    - Outputs the absolute paths of the top-R genomes in the results file, that have >= T/5000 matching
        - Whichever breaks first, T or R
        - How to implement with sequences?
    - Ideally, it should be something like "inspect", where you can look at results and input f.x. '-T 11' and hit Enter.
      Then, the program should print the top 11 results to a file, with full absolute paths separated by a space.
      If you input '-R 4990', the program would print the results having > 4990/5000 in sketch comparison.
      For '-p 0.5', it outputs results having >= 50% of sequences.
      For '-l 2', it outputs results having lost less than or equal to 2 sequences.
      For '-T 11 -R 4990 -p 0.5 -l 2', it would output at most 11 results and, of those, only results satisfying all other constraints.

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