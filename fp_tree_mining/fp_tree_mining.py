"""
Frequent Itemset Mining using FP-tree Algorithm

This program implements the FP-tree growth algorithm to find all frequent itemsets
with support >= 400 in the given dataset. It reads input transaction files and a
vocabulary file, then generates output files containing the frequent itemsets.

Input:
- topic-x.txt: Input transaction files where each line represents a transaction with
  indices of terms.
- vocab.txt: The dictionary mapping term index to term.

Output:
- pattern-x.txt: Output files containing frequent itemsets sorted in descending order
  of support count, with support_count and terms separated by a tab and space, respectively.

Author: Purity Annah Wandondi
University: St. Paul's University
Course: Bachelor of Science in Computer Science.
"""

import os

class TreeNode:
    """Node in the FP-tree."""
    def __init__(self, name, count, parentNode):
        """
        Initialize a TreeNode object.

        Parameters:
        - name: Name of the item
        - count: Count of the item
        - parentNode: Parent node in the tree
        """
        self.name = name  # Name of the item
        self.count = count  # Count of the item
        self.nodeLink = None  # Link similar items
        self.parent = parentNode  # Parent node in the tree
        self.children = {}  # Children of the node

    def inc(self, count):
        """Increase the count of the node."""
        self.count += count

class FPTree:
    """FP-tree for storing itemsets."""
    def __init__(self, transactions, minSup):
        """
        Initialize an FP-tree.

        Parameters:
        - transactions: Dictionary of transactions and their counts
        - minSup: Minimum support threshold
        """
        self.headerTable = {}  # Header table for node links
        self.root = TreeNode('Null', 1, None)  # Root of the tree
        # Count occurrences of items in transactions
        for trans in transactions:
            for item in trans:
                self.headerTable[item] = self.headerTable.get(item, 0) + transactions[trans]
        # Remove items with support less than minSup
        for k in list(self.headerTable.keys()):
            if self.headerTable[k] < minSup: 
                del self.headerTable[k]
        freqItemSet = set(self.headerTable.keys())
        if len(freqItemSet) == 0: return

        # Initialize header table with count and node link
        for k in self.headerTable:
            self.headerTable[k] = [self.headerTable[k], None]
        # Build FP-tree from transactions
        for trans, count in transactions.items():
            localD = {}
            # Filter items based on frequency
            for item in trans:
                if item in freqItemSet:
                    localD[item] = self.headerTable[item][0]
            if len(localD) > 0:
                # Sort items by frequency
                orderedItems = [v[0] for v in sorted(localD.items(), key=lambda p: p[1], reverse=True)]
                # Update tree with ordered items
                self.updateTree(orderedItems, self.root, count)

    def updateTree(self, items, inNode, count):
        """
        Update the tree with a given transaction.
        
        Parameters:
        - items: List of items in the transaction
        - inNode: Current node in the tree
        - count: Count of the transaction
        """
        # Check if the first item in the transaction is already a child node
        if items[0] in inNode.children:
            # If yes, increase the count of the existing node
            inNode.children[items[0]].inc(count)
        else:
            # If not, create a new node and add it as a child node
            inNode.children[items[0]] = TreeNode(items[0], count, inNode)
            # Update the node link in the header table if needed
            if self.headerTable[items[0]][1] == None:
                self.headerTable[items[0]][1] = inNode.children[items[0]]
            else:
                self.updateHeader(self.headerTable[items[0]][1], inNode.children[items[0]])
        
        # Recursively call updateTree for the remaining items in the transaction
        if len(items) > 1:
            self.updateTree(items[1::], inNode.children[items[0]], count)

    def updateHeader(self, nodeToTest, targetNode):
        """
        Update the node link of the header table.
        
        Parameters:
        - nodeToTest: Node to start updating from
        - targetNode: Node to update the node link to
        """
        # Traverse the node link until the end
        while nodeToTest.nodeLink != None:
            nodeToTest = nodeToTest.nodeLink
        # Update the node link to point to the targetNode
        nodeToTest.nodeLink = targetNode

    def mineTree(self, minSup, preFix, freqItemList):
        """
        Mine the FP-tree for frequent itemsets.
        
        Parameters:
        - minSup: Minimum support threshold
        - preFix: Prefix of the current itemset
        - freqItemList: List to store frequent itemsets
        """
        # Sort items in header table by frequency
        bigL = [v[0] for v in sorted(self.headerTable.items(), key=lambda p: p[1][0])]
        # For each item in the header table
        for basePat in bigL:
            # Create a new frequent itemset by adding the current item to the prefix
            newFreqSet = preFix.copy()
            newFreqSet.add(basePat)
            # Append the new frequent itemset to the list
            freqItemList.append(newFreqSet)
            # Find conditional pattern bases for the current item
            condPattBases = self.findPrefixPath(basePat)
            # Create a conditional FP-tree and recursively mine it
            myCondTree = FPTree(condPattBases, minSup)
            if myCondTree.headerTable != {}:
                myCondTree.mineTree(minSup, newFreqSet, freqItemList)

    def findPrefixPath(self, basePat):
        """
        Find conditional pattern base for a given base pattern.
        
        Parameters:
        - basePat: Base pattern for which conditional pattern base needs to be found
        
        Returns:
        - Dictionary containing conditional pattern base along with their counts
        """
        # Get the node corresponding to the base pattern in the header table
        treeNode = self.headerTable[basePat][1]
        condPats = {}  # Initialize an empty dictionary to store conditional pattern base
        while treeNode != None:
            prefixPath = []  # Initialize an empty list to store the prefix path
            # Ascend from the current node to the root, storing node names in the prefix path
            self.ascendTree(treeNode, prefixPath)
            # If the prefix path contains more than one node (excluding the root)
            if len(prefixPath) > 1:
                # Add the prefix path to the conditional pattern base along with its count
                condPats[frozenset(prefixPath[1:])] = treeNode.count
            # Move to the next node in the header table via the node link
            treeNode = treeNode.nodeLink
        return condPats

    def ascendTree(self, leafNode, prefixPath):
        """
        Ascend from a leaf node to the root, storing node names in the prefix path.
        
        Parameters:
        - leafNode: Leaf node from which to ascend
        - prefixPath: List to store the names of nodes in the prefix path
        """
        # If the leaf node has a parent (i.e., not the root)
        if leafNode.parent != None:
            # Append the name of the leaf node to the prefix path
            prefixPath.append(leafNode.name)
            # Recursively ascend to the parent node
            self.ascendTree(leafNode.parent, prefixPath)

def loadDataSet(fileName):
    """
    Load dataset from a file in the input folder.

    Parameters:
    - fileName: Name of the file containing the dataset

    Returns:
    - List of lists representing the dataset
    """
    input_folder = "input"
    file_path = os.path.join(input_folder, fileName)
    dataSet = []  # Initialize an empty list to store the dataset
    with open(file_path, 'r') as file:
        # Read each line from the file, strip whitespace, and split by spaces to get items
        for line in file:
            dataSet.append(line.strip().split())
    return dataSet

def createInitSet(dataSet):
    """
    Create initial set for FP-tree building.

    Parameters:
    - dataSet: List of lists representing the dataset

    Returns:
    - Dictionary where keys are frozensets of transactions and values are their counts
    """
    retDict = {}  # Initialize an empty dictionary to store the initial set
    for trans in dataSet:
        # Convert the transaction list to a frozenset and count occurrences
        retDict[frozenset(trans)] = retDict.get(frozenset(trans), 0) + 1
    return retDict

def loadVocab(fileName):
    """
    Load vocabulary from a file in the input folder.

    Parameters:
    - fileName: Name of the file containing the vocabulary

    Returns:
    - Dictionary mapping term indices (as integers) to terms
    """
    input_folder = "input"
    file_path = os.path.join(input_folder, fileName)
    vocab = {}  # Initialize an empty dictionary to store the vocabulary
    with open(file_path, 'r') as file:
        for line in file:
            # Split each line into index and term, convert index to integer, and add to vocab
            index, term = line.strip().split()
            vocab[int(index)] = term  # Ensure index is converted to an integer
    return vocab

def writePatterns(freqItems, vocab, outputFileName, initSet):
    """
    Write frequent itemsets to file.

    Parameters:
    - freqItems: List of frequent itemsets
    - vocab: Vocabulary dictionary mapping term indices to terms
    - outputFileName: Name of the output file
    - initSet: Initial set for FP-tree building

    Writes the frequent itemsets along with their support counts to the output file.
    """
    patternDict = {}  # Initialize an empty dictionary to store patterns and their counts
    for itemset in freqItems:
        # Calculate support count for each itemset
        supportCount = sum([initSet.get(frozenset(item), 0) for item in initSet if itemset.issubset(item)])
        # Create pattern string by mapping indices to terms from the vocabulary
        pattern = ' '.join(vocab.get(int(item), 'UNKNOWN') for item in itemset if int(item) in vocab)
        patternDict[pattern] = supportCount
    
    # Create the output directory if it doesn't exist
    output_folder = "output"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Write patterns and their counts to the output file
    output_file_path = os.path.join(output_folder, outputFileName)
    with open(output_file_path, 'w') as file:
        for pattern, count in sorted(patternDict.items(), key=lambda x: x[1], reverse=True):
            file.write(f"{count}\t{pattern}\n")

# Main execution function
def main():
    for i in range(1, 5):
        topic_file = f"topic-{i}.txt"
        pattern_file = f"pattern-{i}.txt"
        
        dataSet = loadDataSet(topic_file)
        initSet = createInitSet(dataSet)
        vocab = loadVocab("vocab.txt")
        myFPtree = FPTree(initSet, 400)  # Using 400 as the minimum support
        freqItems = []
        myFPtree.mineTree(400, set([]), freqItems)  # Mining the FP-tree
        writePatterns(freqItems, vocab, pattern_file, initSet)  # Writing the output to pattern-x.txt

if __name__ == "__main__":
    main()
