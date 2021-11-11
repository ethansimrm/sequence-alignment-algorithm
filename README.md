# sequence-alignment-algorithm

This is my solution to a sequence similarity problem formulated as part of the Algorithmic Thinking course, which was organised by Rice University's Department of Computer Science.

Here, I modify and implement the Compute_Global_Alignment_Scores and Compute_Alignment algorithms to calculate both local and global alignments of two sequences using a dynamic programming approach. Rather than calculate all new alignments from scratch and store the maximum, we store previous calculations in an alignment matrix, and use them to facilitate our future alignment scores. The main benefit here is that you won't need to redo the whole calculation if you're just adding dashes!

My work thus far is mirrored in seqalign_alg_implement_mirror.py, and due to the highly specific nature of the modules imported, can only run in CodeSkulptor (a browser-based IDE also created by the Rice University Department of Computer Science) - you can access it at https://py2.codeskulptor.org/#user48_ZImNlQzYEq_57.py.

I then apply these algorithms to the human and fruitfly Eyeless protein sequences, and calculate their similarity via local alignment. To assess if these two local "excerpts" (so to speak) are the PAX domain, I compare these excerpts to a consensus sequence of the PAX domain. I ask whether such similarity could conceivably have arisen by chance, and use null hypothesis testing (human Eyeless vs randomised fruitfly Eyeless) to show that this is highly unlikely. Finally, I use my global alignment algorithm to calculate the edit distance of a string, and use that to compute the set of words within an edit distance of one from the string "humble" and the set of words within an edit distance of two from the string "firefly". 

This second portion is mirrored in similarity_checker_mirror.py, and again, due to the highly specific nature of the modules imported, can only run in CodeSkulptor. It can be accessed at https://py2.codeskulptor.org/#user48_dYdMqaDjOm_6.py.
