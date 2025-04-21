"""
William Schlough
1615902
Homework #1 - PBWT

For HW#1, we read a paper that dealt with the positional Burrows-Wheeler transformation (PBWT) and then completed
the following problem set which directly related back to the paper and previous in-class notebook assignments.
While complex at first - lots of notation and some gibberish that I had to lookup - the coding aspect really helped
me understand what they were trying to explain/convey; the increased speed compared to a naive matching algorithm,
and how we can use transformations to turn a list of strings into a (reverse) prefix array, which we can transform
again and then use said transformation to return to the original list of strings.

Attached are links to some notes for personal use, I always find it helpful to work out what I'm doing on paper.

https://imgur.com/a/37In4DU

https://imgur.com/a/pu1ArVZ
"""


import sys
import numpy

"""The following uses Python to challenge you to create an algorithm for finding
matches between a set of aligned strings. Minimal familiarity with Python is 
necessary, notably list and Numpy array slicing. 
"""

"""Problem 1.

Let X be a list of M binary strings (over the alphabet { 0, 1 }) each of length
N.

For integer 0<=i<=N we define an ith prefix sort as a lexicographic sort
(here 0 precedes 1) of the set of ith prefixes: { x[:i] | x in X }.
Similarly an ith reverse prefix sort is a lexicographic sort of the set of
ith prefixes after each prefix is reversed.

Let A be an Mx(N+1) matrix such that for all 0<=i<M, 0<=j<=N, A[i,j] is the
index in X of the ith string ordered by jth reverse prefix. To break ties
(equal prefixes) the ordering of the strings in X is used.

Complete code for the following function that computes A for a given X.

Here X is a Python list of Python strings.
To represent A we use a 2D Numpy integer array.

Example:

>>> X = getRandomX() #This is in the challenge1UnitTest.py file
>>> X
['110', '000', '001', '010', '100', '001', '100'] #Binary strings, M=7 and N=3
>>> A = constructReversePrefixSortMatrix(X)
>>> A
array([[0, 1, 1, 1],
       [1, 2, 2, 4],
       [2, 3, 5, 6],
       [3, 5, 4, 3],
       [4, 0, 6, 0],
       [5, 4, 3, 2],
       [6, 6, 0, 5]])
>>>

Hint:
Column j (0 < j <= N) of the matrix can be constructed from column j-1 and the
symbol in each sequence at index j-1.

Question 1: In terms of M and N what is the asymptotic cost of your algorithm?

"If we know the sort order at position K ... derive sort order at position (k+1) by a simple process
 looking only at k-th value of each sequence ...
 ... We can therefore calculate the entire set of ordering for all k in a single pass through all the sequences,
in time proportional to NM" (2.1, PBWT)

Based on the reading, the cost of my algorithm is O(NM).
"""

def constructReversePrefixSortMatrix(X):
    #Creates the Mx(N+1) matrix
    A = numpy.empty(shape=[len(X), 1 if len(X) == 0 else len(X[0])+1 ], dtype=int)

    #Code to write - you're free to define extra functions
    #(inline or outside of this function) if you like.

    # A is just a memory block when using numpy.empty, need to fill it in
    # setup columns
    for i in range(len(X)):
        A[i][0] = i

    """# N == number of columns
    N = len(X[0])
    # M == list housing # of rows to iterate over range(0, M-1)
    M = []
    for i in range(len(X)):
        M.append(i)"""
    # follows pseudocode from PBWT paper algorithm #1
    # for each column
    for column in range(len(X[0])):
        # create empty arrays (lists)
        a, b = [], []
        # for row in range 0 to M-1
        for row in range(len(X)):
            # if match
            if X[A[row][column]][column] == '0':
                # append to array a
                a.append(A[row][column])
            else:
                # otherwise, append to array b
                b.append(A[row][column])

        # 3 examples on how to use this style of slicing, last link very good
        # https://stackoverflow.com/questions/509211/understanding-slice-notation
        # https://www.w3schools.com/python/ref_func_range.asp
        # https://stackoverflow.com/questions/33491703/meaning-of-x-x-1-in-python --> 1st example
        # take all rows, the next column will house the concatenation of a + b
        A[ : , column + 1] = a + b # concatenation of a followed by b

        """M = a + b # the concatenation of a followed by b
        # now fill in rows
        for i in range(len(A)):
            A[i][column + 1] = M[i]"""

    return A

"""Problem 2:

Following on from the previous problem, let Y be the MxN matrix such that for
all 0 <= i < M, 0 <= j < N, Y[i,j] = X[A[i,j]][j].

Complete the following to construct Y for X.

Hint: You can either use your solution to constructReversePrefixSortMatrix()
or adapt the code from that algorithm to create Y without using
constructReversePrefixSortMatrix().

Question 2: In terms of M and N what is the asymptotic cost of your algorithm?
Since we use constructReversePrefixSortMatrix(), cost is O(NM), and building Y from
X should yield O(NM) time too.

"""
def constructYFromX(X):
    #Creates the MxN matrix
    Y = numpy.empty(shape=[len(X), 0 if len(X) == 0 else len(X[0]) ], dtype=int)

    #Code to write - you're free to define extra functions
    #(inline or outside of this function) if you like.


    # use constructReversePrefixSortMatrix to start with A (runtime O(NM))
    A = constructReversePrefixSortMatrix(X)

    # note -- Y created as MxN matrix
    # for column in range of all columns in X (which is N-1)
    for column in range(len(X[0])):
        # for row in range M - 1
        for row in range(len(X)):
            # construct Y matrix
            Y[row][column] = X[A[row][column]][column]

    return Y

"""Problem 3.

Y is a transformation of X. Complete the following to construct X from Y,
returning X as a list of strings as defined in problem 1.
Hint: This is the inverse of X to Y, but the code may look very similar.

Question 3a: In terms of M and N what is the asymptotic cost of your algorithm?
Based on problems 1 + 2, asymptotic cost should be O(NM) still. Since these are transformations,
as long as we continue to follow the same procedure as the first problem it should still take
O(NM) time to populate X.


Question 3b: What could you use the transformation of Y for?
Hint: consider the BWT.
We can use Y to go backwards and find the original string / sequence X, similar to the BWT where
we take the last chars of the rotated + sorted string (example used in class that I liked was the cylindrical key-lock)
and then use LF mapping to return the original string.


Question 3c: Can you come up with a more efficient data structure for storing Y?
Personally I can't think of a more efficient data structure for storing Y -- my understanding of numpy is that it is
quite efficient (more so that lists + dictionaries). In the paper's results section + in Tables 1./2., they discuss
the differences between using run-length encoding and gzip to compress the data, and using run-length encoding it compressed by
"more than a factor of a hundred smaller than using gzip," so if there was any way to have a more efficient data
structure it would be to use run-length encoding.

"""
def constructXFromY(Y):
    #Creates the MxN matrix
    X = numpy.empty(shape=[len(Y), 0 if len(Y) == 0 else len(Y[0]) ], dtype=int)

    #Code to write - you're free to define extra functions
    #(inline or outside of this function) if you like.

    # note -- X created as MxN matrix from previous shape of Y
    # Y[i,j] = X[A[i,j]][j]
    # this time, we are creating X from Y
    # from hint, think I can try a variation of problem 1
    # going back to using M as list of indices
    M = []
    for i in range(len(Y)):
        M.append(i)

    # follow similar procedure to problem1
    # for each column
    for column in range(len(Y[0])):
        a, b = [], []
        # for each row in range M-1
        for row in range(len(Y)):
            # here is an addition
            # for each row in range M-1 we want to update current position
            for i in range(len(Y)):
                X[M[i]][column] = Y[i][column]

            # similar to problem 1, we look for match + append to either a or b
            # depending on if the prefix
            if Y[row][column] == 0:
                a.append(M[row])
            else:
                b.append(M[row])

        M = a + b # concatenation of a+b

    # https://www.geeksforgeeks.org/python-map-function/
    # https://docs.aws.amazon.com/lambda/latest/dg/lambda-python.html
    # https://realpython.com/python-lambda/
    # using lambda combined with map(function, iter)
    # --> uses X matrix, and creates a string from each rows entry using ''.join
    return list(map(lambda i : "".join(map(str, i)), X)) #Convert back to a list of strings

"""Problem 4.

Define the common suffix of two strings to be the maximum length suffix shared
by both strings, e.g. for "10110" and "10010" the common suffix is "10" because
both end with "10" but not both "110" or both "010".

Let D be a Mx(N+1) Numpy integer array such that for all 1<=i<M, 1<=j<=N,
D[i,j] is the length of the common suffix between the substrings X[A[i,j]][:j]
and X[A[i-1,j]][:j].

Complete code for the following function that computes D for a given A.

Example:

>>> X = getRandomX()
>>> X
['110', '000', '001', '010', '100', '001', '100']
>>> A = constructReversePrefixSortMatrix(X)
>>> A
array([[0, 1, 1, 1],
       [1, 2, 2, 4],
       [2, 3, 5, 6],
       [3, 5, 4, 3],
       [4, 0, 6, 0],
       [5, 4, 3, 2],
       [6, 6, 0, 5]])
>>> D = constructCommonSuffixMatrix(A, X)
>>> D
array([[0, 0, 0, 0],
       [0, 1, 2, 2],
       [0, 1, 2, 3],
       [0, 1, 1, 1],
       [0, 0, 2, 2],
       [0, 1, 0, 0],
       [0, 1, 1, 3]])

Hints:

As before, column j (0 < j <= N) of the matrix can be constructed from column j-1
and the symbol in each sequence at index j-1.

For an efficient algorithm consider that the length of the common suffix
between X[A[i,j]][:j] and X[A[i-k,j]][:j], for all 0<k<=i is
min(D[i-k+1,j], D[i-k+2,j], ..., D[i,j]).

Question 4: In terms of M and N what is the asymptotic cost of your algorithm?
Based on the reading, even though we are creating a common suffix matrix I'm still are creating
and algorithm with O(NM) asymptotic cost.
"""

def constructCommonSuffixMatrix(A, X):
    D = numpy.zeros(shape=A.shape, dtype=int) #Creates the Mx(N+1) D matrix

    #Code to write - you're free to define extra functions
    #(inline or outside of this function) if you like.

    # note -- D is a Mx(N+1) matrix ARRAY OF ZEROS NOT EMTPY !!!
    # p <- k + 1
    # q <- k + 1
    # for all 0<k<=i -- maybe try setting p + q to 1 at start?
    # for all 1<=j<=N
    for column in range(1, len(A[0])):
        a, b = [], []
        p, q = 0, 0
        # following along with algorithm 2 from paper but we want suffix, no divergent/prefix
        # for i = 0 -> M-1
        for row in range(len(A)):
            # note -- construct from column j-1
            # for 0/1 sort
            # choose min from p, value of D[row][column-1]+1
            p = min(D[row][column - 1] + 1, p)
            # choose min from q, value of D[row][column-1]+1
            q = min(D[row][column - 1] + 1, q)
            # if Y[i][K] == 0
            # note - we take in A as an argument (along with X)
            if X[A[row][column-1]][column-1] == '0':
                a.append(p)
                # have to set p back to 0? NOPE, only when we start at a new column
                # here we want to set p to highest value, allowing min to be either 1 or 0
                # https://www.geeksforgeeks.org/python-infinity/
                p = float('inf')
            else:
                b.append(q)
                # same here, set q to large infinity
                q = float('inf')

        # using same idea from problem 1, view comments above line that looks like this (in P1) for more information
        D[ : , column] = a + b # concatenation of a followed by b

    return D

"""Problem 5.

For a pair of strings X[x], X[y], a long match ending at j is a common substring
of X[x] and X[y] that ends at j (so that X[x][j] != X[y][j] or j == N) that is longer
than a threshold 'minLength'. E.g. for strings "0010100" and "1110111" and length
threshold 2 (or 3) there is a long match "101" ending at 5.

The following algorithm enumerates for all long matches between all substrings of
X, except for simplicity those long matches that are not terminated at
the end of the strings.

Question 5a: What is the asymptotic cost of the algorithm in terms of M, N and the
number of long matches?
From the PBWT paper -- the algorithm in asymptotic cost is O(NM, L) where L is the number of long matches (2.2, 1268)


Question 5b: Can you see any major time efficiencies that could be gained by
refactoring?
While I understand the space efficiency as described by the paper, I think I'm still a little stumped on a major time
efficiency. They describe O(M) space but I don't know if this is the same as reducing asymptotic runtime? If we report
the pairs of subsets a[] and b[] we can run in time O(NM) (2.2, 1268)


Question 5c: Can you see any major space efficiencies that could be gained by
refactoring?
From the PBWT paper -- efficiency can be gained from refactoring and "can be carried out in O(M) space" as long
as we build the A + D while we are running getLongMatches and are "happy to discard previous values of A[i,j] and
D[i,j] for j < k," where i=row, j=column, and k=minlength (2.2, 1268)


Question 5d: Can you imagine alternative algorithms to compute such matches?,
if so, what would be the asymptotic cost and space usage?
From the PBWT paper -- alternative to using suffix / prefix array methods is to build a Hash Table to identify exact
seed matches. While faster due to simpler basic operations, they require more space usage (4, 1271).
Looking back at slides (ESC, 4), using a hash map should give us an asymptotic cost of O(N) where N is the # of seqs,
and the space usage should also be O(N), since a hash map has its space increasing linearly with number of seqs.
"""
def getLongMatches(X, minLength):
    assert minLength > 0
    
    A = constructReversePrefixSortMatrix(X)
    D = constructCommonSuffixMatrix(A, X)
    
    #For each column, in ascending order of column index
    for j in range(1, 0 if len(X) == 0 else len(X[0])):
        #Working arrays used to store indices of strings containing long matches
        #b is an array of strings that have a '0' at position j
        #c is an array of strings that have a '1' at position j
        #When reporting long matches we'll report all pairs of indices in b X c,
        #as these are the long matches that end at j.
        b, c = [], []
        
        #Iterate over the aligned symbols in column j in reverse prefix order
        for i in range(len(X)):
            #For each string in the order check if there is a long match between
            #it and the previous string.
            #If there isn't a long match then this implies that there can
            #be no long matches ending at j between sequences indices in A[:i,j]
            #and sequence indices in A[i:,j], thus we report all long matches
            #found so far and empty the arrays storing long matches.
            if D[i,j] < minLength:
                for x in b:
                    for y in c:
                        #The yield keyword converts the function into a
                        #generator - alternatively we could just to append to
                        #a list and return the list
                        
                        #We return the match as tuple of two sequence
                        #indices (ordered by order in X) and coordinate at which
                        #the match ends
                        yield (x, y, j) if x < y else (y, x, j)
                b, c = [], []
            
            #Partition the sequences by if they have '0' or '1' at position j.
            if X[A[i,j]][j] == '0':
                b.append(A[i,j])
            else:
                c.append(A[i,j])
        
        #Report any leftover long matches for the column
        for x in b:
            for y in c:
                yield (x, y, j) if x < y else (y, x, j)
