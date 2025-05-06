import unittest
import random
import challenge1 as cA

"""The following unitests check the correctness of implementations of 
the functions to be completed.
"""

class TestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.testNo = 100
    
    def testConstructReversePrefixSortMatrix(self):
        """Tests cA.constructReversePrefixSortMatrix using randomly generated 
        examples, checking the output using a simple brute force algorithm.
        """
        for test in range(self.testNo): #Do self.testNo random tests
            X = getRandomX() #A random set of equal length substrings

            A = cA.constructReversePrefixSortMatrix(X)
            
            #Might un-comment the following for debugging
            #print "Test", test, "X", X, "A", A

            #Check has expected number of rows
            self.assertEqual(len(X), A.shape[0]) 
            #Check has expected number of columns
            self.assertEqual(len(X[0])+1 if len(X) > 0 else 1, A.shape[1]) 
            #For each column check the reverse prefix sort
            for j in range(A.shape[1]):
                #Gets a sorted list of reversed prefixes 
                #paired with the index of the parent string in X
                k = sorted([(X[i][:j][::-1], i) for i in range(len(X))]) 
                #The lists of sorted indices should be equal
                self.assertEqual([i[1] for i in k], list(A[:,j]))



    def testConstructYFromX(self):
        """Tests cA.constructXFromY using randomly generated 
        examples.
        """
        for test in range(self.testNo): #Do self.testNo random tests
            X = getRandomX()
            
            #If testConstructReversePrefixSortMatrix is failing then the 
            #inputs to cA.constructCommonSuffixMatrix will be incorrect
            #and this test will likely fail. Hence checking that 
            #cA.constructCommonSuffixMatrix works correctly first is likely
            #to be a good plan.
            A = cA.constructReversePrefixSortMatrix(X)
            Y = cA.constructYFromX(X)
            for j in range(A.shape[1]-1):
                #Check that Y[A[i,j],j] == X[i][j]
                self.assertEqual([int(X[A[i,j]][j]) for i in range(A.shape[0])], list(Y[:,j]))

    
    def testConstructXFromY(self):
        """Tests cA.constructXFromY using randomly generated 
        examples. 
        """
        for test in range(self.testNo): #Do self.testNo random tests
            X = getRandomX()
            Y = cA.constructYFromX(X) #If testConstructYFromX then this will 
            #also likely fail
            X2 = cA.constructXFromY(Y)
    
            self.assertEqual(X, X2)

    def testConstructCommonSuffixMatrix(self):
        """Tests cA.constructCommonSuffixMatrix using randomly generated 
        examples, again checking the output using a simple brute force 
        algorithm.
        """
        for test in range(self.testNo): #Do self.testNo random tests
            X = getRandomX()
            
            #If testConstructReversePrefixSortMatrix is failing then the 
            #inputs to cA.constructCommonSuffixMatrix will be incorrect
            #and this test will likely fail. 
            A = cA.constructReversePrefixSortMatrix(X)
            
            D = cA.constructCommonSuffixMatrix(A, X)
            
            #Might un-comment the following for debugging
            #print "Test", test, "X", X, "D", D
            
            #Check has expected number of rows and columns
            self.assertEqual(A.shape, D.shape) 
            
            #For each column check the length of the common suffixes
            for j in range(A.shape[1]):
                #A sorted list of reversed prefixes
                k = sorted([X[i][:j][::-1] for i in range(len(X))])
                def commonPrefix(x1, x2): 
                    #Computes the common prefix of the two strings x1 and x2
                    i = 0
                    while i < len(x1) and x1[i] == x2[i]:
                        i += 1
                    return x1[:i]
                
                #Convert k into the lengths of the common suffixes of 
                #the strings of interest
                k = ([0] if len(X) > 0 else []) + [len(commonPrefix(k[i], k[i-1])) for i in range(1, len(X))]
                
                #Check the common prefix lengths agree 
                self.assertEqual(k, list(D[:,j]))
    
    def testGetLongMatches(self):
        """Tests cA.getLongMatches using randomly generated
        examples, again checking the output using a simple brute force 
        algorithm.
        """
        for test in range(self.testNo): #Do self.testNo random tests
            X = getRandomX()
            minLength = random.choice(range(1, 10))
            
            #Get the matches from cA.getAllLongMatches putting them into a 
            #list/set
            matchesList = list(cA.getLongMatches(X, minLength))
            matches = set(matchesList) 
            
            #Check we don't have any duplicates in the output
            self.assertEqual(len(matches), len(matchesList)) 
            
            #Very simple brute force algorithm to compute long matches
            matches2 = set()
            for i in range(len(X)):
                x1 = X[i]
                for k in range(i+1, len(X)):
                    x2 = X[k]
                    matchLength = 0
                    for j in range(len(x1)):
                        if x1[j] == x2[j]:
                            matchLength += 1
                        else:
                            if matchLength >= minLength:
                                matches2.add((i, k, j))
                            matchLength = 0
            
            #Check the sets are equivalent
            self.assertEqual(matches2, matches)

def getRandomX():
    """Creates a small random set of binary strings all of the same length.
    """
    seqNo = random.choice(range(1, 10)) 
    seqLength = random.choice(range(1, 10))
    randChar = lambda : "0" if random.random() > 0.5 else "1"
    randSeq = lambda : "".join([ randChar() for i in range(seqLength)])
    ###had to rewrite this function because return map() doesn't work in python 3 (returns a map object and not your actual set of strings)
    return [randSeq() for _ in range(seqNo)]

if __name__ == '__main__':
    unittest.main()