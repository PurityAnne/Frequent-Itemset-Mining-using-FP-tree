# Frequent-Itemset-Mining-using-FP-tree

This program implements the FP-tree growth algorithm to find all frequent itemsets with support >= 400 in the given dataset. It reads input transaction files and a vocabulary file, then generates output files containing the frequent itemsets.

Input:
- topic-x.txt: Input transaction files where each line represents a transaction with indices of terms.
- vocab.txt: The dictionary mapping term index to term.

Output:
- pattern-x.txt: Output files containing frequent itemsets sorted in descending order   of support count, with support_count and terms separated by a tab and space, respectively.
