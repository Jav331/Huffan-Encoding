/*
 * This program implements a compression and decompression
 * via a huffman tree algorithm.
 */
#include "bits.h"
#include "treenode.h"
#include "huffman.h"
#include "map.h"
#include "vector.h"
#include "priorityqueue.h"
#include "strlib.h"
#include "SimpleTest.h"

using namespace std;

/**
 * Given a Queue<Bit> containing the compressed message bits and the original
 * encoding tree, decode the bits back to the original message text.
 *
 * You can assume that tree is a well-formed non-empty encoding tree and the
 * messageBits queue contains a valid sequence of encoded bits for that tree.
 *
 * Your implementation may change the messageBits queue however you like. There
 * are no requirements about its contents after this function returns.
 * The encoding tree should be unchanged.
 *
 */
string decodeText(EncodingTreeNode* tree, Queue<Bit>& messageBits) {


    string result;
    EncodingTreeNode* treePath = tree;

    while (!messageBits.isEmpty()) {
        if (messageBits.peek() == 0) {
            messageBits.dequeue();
            treePath = treePath->zero;
        } else if (messageBits.peek() == 1) {

            messageBits.dequeue();
            treePath = treePath->one;
        }

        if (treePath->isLeaf()) {
            result += treePath->getChar();
            // go back to top of tree
            treePath = tree;
        }
    }



    return result;
}

/**
 * Reconstruct an encoding tree from flattened form Queue<Bit> and Queue<char>.
 *
 * You can assume that the queues are well-formed and the shape/leaves sequences
 * form a valid encoding tree.
 *
 * Your implementation may change the treeShape and treeLeaves queues however you like.
 * There are no requirements about their contents after this function returns.
 *
 */
EncodingTreeNode* unflattenTree(Queue<Bit>& treeShape, Queue<char>& treeLeaves) {

    if (treeShape.isEmpty()) {
        return nullptr;
    }


    // if 0, this is our leaf
    if (treeShape.peek() == 0) {
        treeShape.dequeue();
        char leafNode = treeLeaves.dequeue();
        EncodingTreeNode *result = new EncodingTreeNode(leafNode);
        return result;

    // treeShape.peek() == 1
    } else { // keep traversing
        treeShape.dequeue();
        EncodingTreeNode *leftSub = unflattenTree(treeShape, treeLeaves);
        EncodingTreeNode *rightSub = unflattenTree(treeShape, treeLeaves);
        EncodingTreeNode *result = new EncodingTreeNode(leftSub, rightSub);
        return result;
    }

}

/*
 * Decompress the given EncodedData and return the original text.
 *
 * You can assume the data argument is well-formed and was created by a correct
 * implementation of compress.
 *
 * Your implementation may change the EncodedData however you like. There
 * are no requirements about its contents after this function returns.
 *
 */
string decompress(EncodedData& data) {

    string result;

    EncodingTreeNode *tree = unflattenTree(data.treeShape, data.treeLeaves);

    result = decodeText(tree, data.messageBits);

    // deallocate tree when done
    deallocateTree(tree);

    return result;
}
/*
* Helper function that builds a frequency map
*/
void getFrequencies(string text, Map<char, int> &map) {
    for (char ch: text) {
        map[ch]++;
    }
}

/**
 * Constructs an optimal Huffman coding tree for the given text, using
 * the Huffman algorithm.
 *
 * Reports an error if the input text does not contain at least
 * two distinct characters.
 *
 * It can be helpful for you to establish/document the expected behavior of
 * choices re: tie-breaking, choice of which subtree on which side.
 * These choices do not affect correctness or optimality of resulting tree
 * but knowing which is used allows you to construct test cases that depend
 * on that behavior. Our provided test cases expect:
 *  -- our pqueue dequeues elems of equal priority in FIFO order
 *  -- when build new interior node, first elem dequeued placed as zero subtree
 *     and second elem dequeued placed as one subtree
 *
 *   referenced Sets and Maps lecture to build
 *   getFrequencies function
 */
EncodingTreeNode* buildHuffmanTree(string text) {

    // get frequency of char in text
    Map<char, int> freqChar;
    getFrequencies(text, freqChar);

    if (freqChar.size() < 2) {
        error("Input text must contain at least two distinct characters in order to build a Huffman encoding.");
    }

    PriorityQueue<EncodingTreeNode*> huffmanPQ;

    for (char ch: freqChar) {
        huffmanPQ.enqueue(new EncodingTreeNode(ch), freqChar[ch]);
    }

    while (huffmanPQ.size() > 1) {

        int priorityTemp1 = huffmanPQ.peekPriority();
        EncodingTreeNode *temp1 = huffmanPQ.dequeue();
        int priorityTemp2 = huffmanPQ.peekPriority();
        EncodingTreeNode *temp2 = huffmanPQ.dequeue();


        int priorityFreq = priorityTemp1 + priorityTemp2;

        EncodingTreeNode *connectNode = new EncodingTreeNode(temp1, temp2);

        huffmanPQ.enqueue(connectNode, priorityFreq);

    }

    EncodingTreeNode *huffmanTree = huffmanPQ.dequeue();


    return huffmanTree;
}


void buildTableHelper (EncodingTreeNode* tree, Map<char, Queue<Bit>> &table, Queue<Bit> &path) {

    if (tree->isLeaf()) {
        table[tree->getChar()] = path;
    } else if (tree->zero != nullptr){

        // choose -> explore -> unchose

        // left tree
        Queue<Bit> left = path;
        left.enqueue(0);
        buildTableHelper(tree->zero, table, left);

        if (tree->one != nullptr) {

            // right tree
            Queue<Bit> right = path;
            right.enqueue(1);
            buildTableHelper(tree->one, table, right);

        }
    }

}

Map<char, Queue<Bit>> buildTable(EncodingTreeNode* tree) {

    Map<char, Queue<Bit>> table;
    Queue<Bit> path;

    buildTableHelper (tree, table, path);

    return table;
}


/**
 * Given a string and an encoding tree, encode the text using the tree
 * and return a Queue<Bit> of the encoded bit sequence.
 *
 * You can assume tree is a well-formed encoding tree and contains
 * an encoding for every character in the text.
 *
 */
Queue<Bit> encodeText(EncodingTreeNode* tree, string text) {

    Queue<Bit> resultBits;

    Map<char, Queue<Bit>> table = buildTable(tree);

    for (char ch: text) {
        Queue<Bit> seqForChar = table[ch];
        while(!seqForChar.isEmpty()) {
            resultBits.enqueue(seqForChar.dequeue());
        }
    }

    return resultBits;
}

/**
 * Flatten the given tree into a Queue<Bit> and Queue<char> in the manner
 * specified in the assignment writeup.
 *
 * You can assume the input queues are empty on entry to this function.
 *
 * You can assume tree is a well-formed encoding tree.
 *
 */
void flattenTree(EncodingTreeNode* tree, Queue<Bit>& treeShape, Queue<char>& treeLeaves) {

    if (tree == nullptr) {
        return;
    }


    if (tree->isLeaf()) {
        treeShape.enqueue(0);
        treeLeaves.enqueue(tree->getChar());
    } else {
        treeShape.enqueue(1);
        flattenTree(tree->zero, treeShape, treeLeaves);
        flattenTree(tree->one, treeShape, treeLeaves);
    }}

/**
 * Compress the message text using Huffman coding, producing as output
 * an EncodedData containing the encoded message bits and flattened
 * encoding tree.
 *
 * Reports an error if the message text does not contain at least
 * two distinct characters.
 *
 */
EncodedData compress(string messageText) {


    EncodedData data;
    Queue<Bit>  treeShape;
    Queue<char> treeLeaves;
    Queue<Bit> messageBits;

    EncodingTreeNode* huffmanTree = buildHuffmanTree(messageText);
    messageBits = encodeText(huffmanTree, messageText);
    flattenTree(huffmanTree, treeShape, treeLeaves);

    data.treeLeaves = treeLeaves;
    data.treeShape = treeShape;
    data.messageBits = messageBits;

    deallocateTree(huffmanTree);

    return data;
}

/* * * * * * Testing Helper Functions Below This Point * * * * * */

EncodingTreeNode* createExampleTree() {
    /* Example encoding tree used in multiple test cases:
     *                *
     *              /   \
     *             T     *
     *                  / \
     *                 *   E
     *                / \
     *               R   S
     */

    EncodingTreeNode *t = new EncodingTreeNode('T');
    EncodingTreeNode *r = new EncodingTreeNode('R');
    EncodingTreeNode *e = new EncodingTreeNode('E');
    EncodingTreeNode *s = new EncodingTreeNode('S');

    // levels denote the level of the BST
    // worked from bottom up
    EncodingTreeNode *level3 = new EncodingTreeNode(r, s);
    EncodingTreeNode *level2 = new EncodingTreeNode(level3, e);
    EncodingTreeNode *level1 = new EncodingTreeNode(t, level2);

    EncodingTreeNode *result = level1;



    return result;
}

void deallocateTree(EncodingTreeNode* t) {
    // preformed recursively
    if (t == nullptr) {
        return;
    }

    deallocateTree(t->one);
    deallocateTree(t->zero);
    delete t;
}


bool areEqual(EncodingTreeNode* a, EncodingTreeNode* b) {


    if (a == nullptr && b == nullptr) {
        return true;
    } else if (a == nullptr || b == nullptr){
        return false;
    }

    if (a->isLeaf() && b->isLeaf()) {
        return (a->getChar() == b->getChar());
    }



    bool leftSubTree = areEqual(a->zero, b->zero);
    bool rightSubTree = areEqual(a->one, b->one);

    return leftSubTree && rightSubTree;
}

/* * * * * * Test Cases Below This Point * * * * * */

/* TODO: Write your own student tests. */


STUDENT_TEST("areEqual on empty tree") {
    EncodingTreeNode* a = nullptr;
    EncodingTreeNode* b = nullptr;

    EXPECT(areEqual(a, b));
}

STUDENT_TEST("areEqual on simple tree and create tree") {

    /* Example encoding tree used in multiple test cases:
 *
 *                *
 *              /   \
 *             A     B
 *
 */


    EncodingTreeNode *a = new EncodingTreeNode('A');
    EncodingTreeNode *b = new EncodingTreeNode('B');

    EncodingTreeNode *simpleTree = new EncodingTreeNode(a, b);
    EncodingTreeNode *tree = createExampleTree();

    EXPECT(!areEqual(simpleTree, tree));
    deallocateTree(simpleTree);
    deallocateTree(tree);

}

STUDENT_TEST("areEqual on simple non-empty tree") {

    /* Example encoding tree used in multiple test cases:
 *
 *                *
 *              /   \
 *             A     B
 *
 */


    EncodingTreeNode *a = new EncodingTreeNode('A');
    EncodingTreeNode *b = new EncodingTreeNode('B');

    EncodingTreeNode *simpleTree = new EncodingTreeNode(a, b);
    EncodingTreeNode *emptyTree = nullptr;

    EXPECT(!areEqual(simpleTree, emptyTree));
    deallocateTree(simpleTree);
}


STUDENT_TEST("areEqual on simple two different trees") {


    EncodingTreeNode *a = new EncodingTreeNode('A');
    EncodingTreeNode *b = new EncodingTreeNode('B');
    EncodingTreeNode *simpleTree = new EncodingTreeNode(a, b);


    EncodingTreeNode *c = new EncodingTreeNode('C');
    EncodingTreeNode *d = new EncodingTreeNode('D');
    EncodingTreeNode *e = new EncodingTreeNode('E');
    EncodingTreeNode *f = new EncodingTreeNode('F');
    EncodingTreeNode *tree1 = new EncodingTreeNode(c, d);
    EncodingTreeNode *tree2 = new EncodingTreeNode(tree1, e);
    EncodingTreeNode *tree3 = new EncodingTreeNode(tree2, f);


    EXPECT(!areEqual(simpleTree, tree3));
    deallocateTree(simpleTree);
    deallocateTree(tree3);

}

STUDENT_TEST("areEqual on two same trees") {
    EncodingTreeNode* simpleTree1 = createExampleTree();
    EncodingTreeNode* simpleTree2 = createExampleTree();

    EXPECT(areEqual(simpleTree1, simpleTree2));

    deallocateTree(simpleTree1);
    deallocateTree(simpleTree2);
}

STUDENT_TEST("deallocateTree") {
    EncodingTreeNode* tree = createExampleTree();
    deallocateTree(tree);

}






/* * * * * Provided Tests Below This Point * * * * */



PROVIDED_TEST("decodeText, small fixed inputs, example encoding tree") {
    EncodingTreeNode* tree = createExampleTree(); // see diagram above
    EXPECT(tree != nullptr);

    Queue<Bit> messageBits = {1,1}; // E
    EXPECT_EQUAL(decodeText(tree, messageBits), "E");

    messageBits = {1,0,1,1,1,0}; // SET
    EXPECT_EQUAL(decodeText(tree, messageBits), "SET");

    messageBits = {1,0,1,0,1,0,0,1,1,1,1,0,1,0,1}; // STREETS
    EXPECT_EQUAL(decodeText(tree, messageBits), "STREETS");

    deallocateTree(tree);
}

PROVIDED_TEST("unflattenTree, example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above
    Queue<Bit>  treeShape  = {1,0,1,1,0,0,0};
    Queue<char> treeLeaves = {'T','R','S','E'};
    EncodingTreeNode* tree = unflattenTree(treeShape, treeLeaves);

    EXPECT(areEqual(tree, reference));

    deallocateTree(tree);
    deallocateTree(reference);
}

PROVIDED_TEST("decompress, small fixed input, example encoding tree") {
    EncodedData data = {
        {1,0,1,1,0,0,0},            // treeShape
        {'T','R','S','E'},          // treeLeaves
        {0,1,0,0,1,1,1,0,1,1,0,1}   // messageBits
    };

    EXPECT_EQUAL(decompress(data), "TRESS");
}

PROVIDED_TEST("buildHuffmanTree, small fixed input, example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above
    EncodingTreeNode* tree = buildHuffmanTree("STREETTEST");
    EXPECT(areEqual(tree, reference));

    deallocateTree(reference);
    deallocateTree(tree);
}

PROVIDED_TEST("encodeText, small fixed inputs, example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above

    Queue<Bit> messageBits = {1,1}; // E
    EXPECT_EQUAL(encodeText(reference, "E"), messageBits);

    messageBits = {1,0,1,1,1,0};    // SET
    EXPECT_EQUAL(encodeText(reference, "SET"), messageBits);

    messageBits = {1,0,1,0,1,0,0,1,1,1,1,0,1,0,1}; // STREETS
    EXPECT_EQUAL(encodeText(reference, "STREETS"), messageBits);

    deallocateTree(reference);
}

PROVIDED_TEST("flattenTree, example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above
    Queue<Bit>  expectedShape  = {1,0,1,1,0,0,0};
    Queue<char> expectedLeaves = {'T','R','S','E'};

    Queue<Bit>  treeShape;
    Queue<char> treeLeaves;
    flattenTree(reference, treeShape, treeLeaves);

    EXPECT_EQUAL(treeShape,  expectedShape);
    EXPECT_EQUAL(treeLeaves, expectedLeaves);

    deallocateTree(reference);
}

PROVIDED_TEST("compress, small fixed input, example encoding tree") {
    EncodedData data = compress("STREETTEST");
    Queue<Bit>  treeShape   = {1,0,1,1,0,0,0};
    Queue<char> treeChars   = {'T','R','S','E'};
    Queue<Bit>  messageBits = {1,0,1,0,1,0,0,1,1,1,1,0,0,1,1,1,0,1,0};

    EXPECT_EQUAL(data.treeShape, treeShape);
    EXPECT_EQUAL(data.treeLeaves, treeChars);
    EXPECT_EQUAL(data.messageBits, messageBits);
}

PROVIDED_TEST("Test end-to-end compress -> decompress, small fixed inputs") {
    Vector<string> inputs = {
        "HAPPY HIP HOP",
        "Nana Nana Nana Nana Nana Nana Nana Nana Batman",
        "Research is formalized curiosity. It is poking and prying with a purpose. – Zora Neale Hurston",
    };

    for (string input: inputs) {
        EncodedData data = compress(input);
        string output = decompress(data);

        EXPECT_EQUAL(input, output);
    }
}
