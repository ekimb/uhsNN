import math
class Graph:
    def __init__(self, k):
        self.edgeNum = int(math.pow(4, k))
        self.k = k
        self.edgeArray = [None] * self.edgeNum
        self.vertexExp = int(math.pow(4, self.k - 1))

        for i in range(self.edgeNum):
            self.edgeArray[i] = 1

        self.edgeCount = self.edgeNum
        self.vertexCount = self.edgeNum / 4
        self.alphabetMap = {}
        self.alphabet = "ACGT"

        for i in range(4):
            self.alphabetMap[self.alphabet[i]] = i

    def getIndex(self, label):
        index = 0
        for j in range(self.k):
            index = index + self.alphabetMap[label[j]] * int(math.pow(4, self.k - j - 1))
        return index
        
    def getLabel(self, index):
        if index < 0:
            return ''
        finalString = ""
        for j in range(self.k):
            finalString = self.alphabet[(index % 4)] + finalString
            index = int(index / 4)
        return finalString
    
    def topologicalSort(self):
        for i in range(self.vertexExp):
            self.used[i] = False
            self.finished[i] = False
            self.topoSort[i] = 0
        index = 0
        for i in range(self.vertexExp):
            if self.used[i] == False:
                index = self.depthFirstSearch(index, i)
                if (index == -1):
                    self.topoSort = None
                    break
    
    def maxLength(self):
        depth = [0] * self.vertexExp
        maxDepth = -1
        for i in range(self.vertexExp):
            depth[i] = 0
            maxVertDepth = -1
            for j in range(4):
                edgeIndex = self.topoSort[i] + j * self.vertexExp
                vertexIndex = int(edgeIndex / 4)
                if ((depth[vertexIndex] > maxVertDepth) and (self.edgeArray[edgeIndex] == 1)):
                    maxVertDepth = depth[vertexIndex]
            depth[self.topoSort[i]] = maxVertDepth + 1
            if (depth[self.topoSort[i]] > maxDepth):
                 maxDepth = depth[self.topoSort[i]]
        return maxDepth
    
    def depthFirstSearch(self, index, u):
        self.used[u] = True
        cycle = False
        for v in self.getAdjacent(u):
            if self.used[v] == True and self.finished[v] == False:
                cycle = True
            if self.used[v] == False:
                index = self.depthFirstSearch(index, v)
                cycle = cycle or (index == -1)
        self.finished[u] = True
        self.topoSort[index] = u
        if cycle:
            return -1
        else:
            return index + 1

    def getAdjacent(self, v):
        count = 0
        adjVertex = [0] * 4
        for i in range(4):
            index = v + i * self.vertexExp
            if self.edgeArray[index] == 1:
                adjVertex[count] = int(index / 4)
                count = count + 1
        rc = []
        for i in range(count):
            rc.append(adjVertex[i])
        return rc

    def removeEdge(self, i):
        self.edgeCount = self.edgeCount - 1
        self.edgeArray[i] = 0
    
    def updateNumbers(self, L):
        maxHittingNum = 0
        imaxHittingNum = -1
        for i in range(self.edgeNum):
            hittingNum = 0
            for j in range(((1 - self.edgeArray[i]) * L), L):
                hittingNum = hittingNum + self.Fcurr[i % self.vertexExp] * self.D[(L-j-1)][int(i / 4)]
                self.hittingNumArray[i] = hittingNum
            if (hittingNum > maxHittingNum):
                maxHittingNum = hittingNum
                imaxHittingNum = i
        return imaxHittingNum
            
    def HittingPipeline(self, L):
        imaxHittingNum = -1
        count = 0
        l = L - self.k + 1
        self.hittingNumArray = [0] * self.edgeNum
        self.used = [0] * self.vertexExp
        self.finished = [0] * self.vertexExp
        self.topoSort = [0] * self.vertexExp
        self.D = []
        for i in range(l+1):
            Dpool = [0] * ((l+1)* self.vertexExp)
            self.D.append(Dpool)
        self.Fprev = [0] * self.vertexExp
        self.Fcurr = [0] * self.vertexExp
        while (self.calculatePaths(l) == True):
            imaxHittingNum = self.updateNumbers(l)
            if imaxHittingNum < 0:
                break
            print(imaxHittingNum)
            self.removeEdge(imaxHittingNum)
            label = self.getLabel(imaxHittingNum)
            print(label)
            count = count + 1
        print('Additional set size: ' + str(count))
        self.topologicalSort()
        print("Length of longest remaining path: " + str(self.maxLength())  + "\n")
        return count

    def calculatePaths(self, L):
        curr = 1
        vertexExp2 = self.vertexExp * 2
        vertexExp3 = self.vertexExp * 3
        vertexExpMask = self.vertexExp - 1
        for i in range(self.vertexExp):
             self.D[0][i] = 1
             self.Fprev[i] = 1
            
        for j in range(1, L+1):
            for i in range(self.vertexExp):
                self.D[j][i] = self.edgeArray[i]*self.D[j-1][(i >> 2)] + self.edgeArray[i + self.vertexExp]*self.D[j-1][((i + self.vertexExp) >> 2)] + self.edgeArray[i + vertexExp2]*self.D[j-1][((i + vertexExp2) >> 2)] + self.edgeArray[i + vertexExp3]*self.D[j-1][((i + vertexExp3) >> 2)]
        for i in range(self.edgeNum):
            self.hittingNumArray[i] = 0
        while (curr <= L):
            for i in range(self.vertexExp):
                index = (i * 4)
                self.Fcurr[i] = (self.edgeArray[index]*self.Fprev[index & vertexExpMask] + self.edgeArray[index + 1]*self.Fprev[(index + 1) & vertexExpMask] + self.edgeArray[index + 2]*self.Fprev[(index + 2) & vertexExpMask] + self.edgeArray[index + 3]*self.Fprev[(index + 3) & vertexExpMask])
            for i in range(self.edgeNum):
                self.hittingNumArray[i] = self.hittingNumArray[i] + (self.Fprev[i % self.vertexExp] * self.D[(L-curr)][int(i / 4)])
                if (self.edgeArray[i] == 0):
                    self.hittingNumArray[i] = 0
            for i in range(self.vertexExp):
                self.Fprev[i] = self.Fcurr[i]
            curr = curr + 1
        return True

    def removeDecycling(self, decycPath):
        with open(decycPath) as f:
            for line in f:
                kmer = line.rstrip()
                self.removeEdge(self.getIndex(kmer))
        print("Decycling set done!")


