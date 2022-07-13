## Tasks

1. Build the algorithm to iterate over sequences from .mince file
    - Needs to be modular, i.e. be able to take in multiple .mince files and not work linearly through them
    - Start off with model in Python
    - Then translate to C++
    - KIND OF DONE?
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