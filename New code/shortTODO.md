## Tasks

### Actually next up:
    - Launched ./trialSequences at around 18:00, Thursday
    - It's an updated version of the old ./sequences, which doesn't get bogged down in hashing vectorForms and the like
    - Removed fasta_mapping from the code, so now the fasta sequence IDs are paired with the vectors in Uniques
    - As soon as that finishes, I run ./hashSequences and create my first fully fledged MinCE library.

    - Why is cliques/118331_0.clique cycling?

    ./NewCode/local/bin/trialSequences -k 31 -N 10 cliques/41620_00.clique
    ./NewCode/local/bin/trialSequences -k 31 -N 10 cliques/46178_02.clique
    ./NewCode/local/bin/trialSequences -k 31 -N 10 cliques/163920_0.clique

#### Where I'm at:
    - I've modified most files within MinCE in some way or another. I've set UFClusters to work on clustering with t=4997.
    - When UFClusters finishes, I should run bottomUpCliques on the resulting clusters.
        - While UFClusters is running, I should take the time to re-examine the architecture of bottomUpCliques
    - When cliques are ready, run Sequence on those clusters. (make trial runs w.r.t. new code for FastX parser)
    - When sequences are ready, run makeSeqHashmap on all_sequences.txt
    - Then try Mince2.cpp, I guess...

#### Problems persist:
    - First of all, bottomUpCliques is likely insufficient. The joining process is too arbitrary.
    - Second, the greedy algorithm in Sequences.cpp is not water-proof, as it may not distinguish between all members.
        - To demonstrate: 
            - Let's say member 3 is on par with every other w.r.t. connections matrix.
            - It is not chosen as 'focus_member', instead we choose 1 and choose vector form [1,0,1,0,0].
            - Next up we choose member 5 as focus_member and choose vector [1,0,1,0,1].
            - Now we've chosen two vectors for member 3, neither of which discerns it from member 1.
            - For larger cliques, the possibility of this number reaching 10 is not out of the question.

#### Next up:
    - I need to create a new program which takes in all .sequence files + old index file and:
        - creates a .hashmap file called seq.hashmap which maps:
            - hash(sequence) -> IDs of genomes having that sequence
        - creates a new .index file called MinCE_to_NCBI.index with format:
            - { MinCE_index         NCBI_index          Number_of_sequences }
            - NCBI index should be trimmed so as only to contain NCBI needed info
            - Number_of_sequences should be -1 for genomes not in a clique
    - Rename hash_locator to sketch.hashmap
    - At final stage, only files needed should be:
        - MinCE_to_NCBI.index
        - sketch.hashmap
        - seq.hashmap
        - ./mince executable

0.  Implement Lior's ideas:
    #### a. Create a reverse dictionary, name_to_ind, mapping { name_of_sketch -> sketch_ind }
    #### b. Create a dictionary, seq_locator, for all sequences in the set, mapping { hash(sequence) -> sketch_ind }
    #### c. Implement MinCE to incorporate seq_locator while running triangulate.
        - Iterate over input fastq file, reading every kmer and turning it into a hash value.
        - Check whether hash value is less than MAX_HASH, in which case log points in res.sketch through hash_locator.
        - Next, send the hash value into seq_locator. If it exists, log points in res.seq through seq_locator.
        - Load index and fetch names of genomes in results, fill into res.name.
        - Sort results by res.seq first and then res.sketch, in event of tie.
        - Print results to terminal.
        - Log result to file, easily read by f.x. kallisto.
    #### d. I think there's room to change the architecture of hash_locator, seq_locator and index here. They seem too many and large.
    #### e. Finish work on sendToCliques.cpp and incorporation with MinCE, if I want an easier time testing stuff?
        - Just finish work on sendToCliques.cpp and then make a copy of Mince.cpp (ModularMince.cpp).
        - Modify the copy to run according to original idea by including sendToCliques.cpp
        - Make a new bin ./modular_mince, use that as a benchmark for the debug of second stage

### ATH: Do I need to change something about batch files as input, now that I've trimmed the fat off everything in Mince2.cpp?

1. Build the algorithm to iterate over sequences from .mince file
    - Needs to be modular, i.e. be able to take in multiple .mince files and not work linearly through them
    - Start off with model in Python
    - Then translate to C++
    - KIND OF DONE?
        - Question: For cliques with no discerning k-mers, do I just output the members of that clique at result stage with 'Likely'?
2. Finish work on Mince.cpp, to make it tailored to big fastq files
    - Multiple results is already pretty much implemented
    - Megasketch needs some testing to determine whether its implementation is concrete
    - I think the biggest question is how we return the results to users
        - In particular, how we deliver singleton results mixed in with results from .minced files
3. Translate bottomUpCliques.py into C++, if it's deemed worth my while
    - It could be attached to UFClusters and compiled into one step, which is a plus.
4. Design a pipeline in SnakeMake to generate a MinCE library, given some parameters
    - Iterate over given parameters to try to find best combination 
    - Or iterate over them to find bugs and get a better understanding of precision
5. Find old articles (maybe 7-8 years old) about metagenomics research on bacteria
    - It might be good to run some workflows on the same datasets as them, to show MinCE's capabilities
    - Kristján or Lior might know more about where to look


## Bug report

./distinguish halts on cliques/224345_0.clique due to Segmentation fault. 
It's a clique with 2 members and does not meet N = 10 for either of them.

Fixed a whole bunch of things, it should work now.