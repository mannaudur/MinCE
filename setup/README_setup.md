# MinCE - setting up your own databank
MinCE - Minhash with Cluster Extensions


This is an instructions manual on how to setup your own databank with MinCE. As of now, this is a fairly rough draft but the basic steps are all there. You're also free to contact us with any questions!

Load your genome fasta files into a directory of your choosing, such as 'toyset/'.

Create another directory for sketches, such as 'toyset/toysketches/'.

Edit massmash.sh to match your directory names and run it. This will take some time.

The run produces output based default values for smash, for example k-mer size k=31 and sketch size s=1000. This can be changed by editing massmash.sh based on the smash syntax.

Next transfer all names in your sketch directory to a .txt file, for example by piping. It is recommended using find instead of ls, as it's buffer is larger: 'find toyset/toysketches/ -name "*.sketch" > all.txt'.

Now create a folder called 'atoms/' and run 'smash/bin/atom all.txt' from its parent directory. This will produce outputs with a default threshold value of > 995/1000 for atom joining, but this can be changed with -l. Note that the paths in 'all.txt' will need to have a proper relevant path to the sketches, from the directory you run 'atom' from.

This will create an extensionless file in 'atoms/' for every atom discovered in your dataset, as well as two reference files; 'hashlocator' and 'indices'. We want to add the extension .txt to the extensionless atom files, svo run rename.sh in the parent directory of 'atoms/'. As we said, it's a rough draft.

The file 'hashlocator' contains every single hash value in your dataset and maps it to every sketch containing it. The sketches are referenced by number, not name, and 'indices' maps the number as an index to genome name and its corresponding 
atom. Atom names are numbers, based on the index value of the parent node from the creation of atoms using union-find. If the sketch has no corresponding atom, its atom value in 'indices' will be NULL.

Next you will need to install Bifrost (https://github.com/pmelsted/bifrost) and add it to your system's $PATH. If you do not have access to your system's $PATH, you can choose to download the binary for Bifrost into another folder. This step is explained on the Bifrost Git, but you will also need to edit 'get_feats.sh' to provide the new path to the binary.

Finally, you copy 'extract_features.py' and 'deBruijn_extractor.py' to the parent directory of 'atom/'. In that same directory, you create another directory called 'features/'.

Now, if you're not worried about overloading your system, run 'get_feats.sh'. This will take a long time, but you can just let it 
simmer.

As an alternative, this can be broken up into seperate processes. We provide two bash scripts to break this into something more managable; 'just_get_bitmats.sh' takes a .txt file of paths and translates them into bitmatrix and fasta files; and 'just_get_feats.sh' takes a .txt with paths to those bitmatrix/fasta files and translates them into .json files containing relevant 
features.

By controlling the .txt files, you can schedule chunks to be translated at each time. Just note that 'just_get_feats.sh' can't run on a bitmatrix/fasta file, unless 'just_get_bitmats.sh' has finished creating it. Generally, 'just_get_bitmats.sh' also takes longer than 'just_get_feats.sh'.

As a last adjustment, go ahead and rename 'indices' to 'index_mapping.txt'. Then run 'get_revDicts.py hashlocator' in the directory containing 'hashlocator'. This should create a folder called 'reverseDicts/' and and 90 subdictionaries therein - one for every pair of starting numbers. For example, 'revDict10.txt' will contain all hash values starting with 10... and 'revDict99.txt' will contain all hash values starting with 99....

Lastly, if you intend on using the Megasketch function, you will need to inspect the value in 'MAX_HASH.txt'. This is your MAX_HASH value, surprisingly. You will need to insert that into the file  'smash/src/Minhash.cpp'. Find the line 

'const uint64_t MAX_HASH = 9999999776999205UL;'

and replace it with

'const uint64_t MAX_HASH = \<MAX_HASH\>UL;'

and then run 'make' again in the 'smash' directory. Note that this value is dataset-specific, so you'll need to adjust to whichever 
dataset you're using, which is why you'll need to update the one we assigned to the original 285K dataset.

That's it, you're done! Just make sure the folders follow the original arrangement and you should be able to run any query on your 
own, personal dataset.

