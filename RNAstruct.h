#ifndef RNASTRUCTURE_H
#define RNASTRUCTURE_H


#include"set.h"
#include"vector.h"
#include"map.h"

using namespace std;
class RNAStructure
{
public:
    RNAStructure(ifstream& input);
    int nextIndex(int index, int prior);


private:
    set<int> m_loopPair;// all sequences frequently found in loops
    vector< char>  m_sequence;//all nucleotides in the sequence, position same as the index in the vector
    map<char, char> m_antisense;// stores values and their antisense for easy location of complements
    map<int, int> m_basePairs;// stores indeces of paired sequences
    size_t m_length;//length of the sequence
    void initiate (ifstream& input, set<int>& intePositions);
    void generate(ifstream& input);
    int findSymmetry(int index, int start, int end);

};



#endif // RNASTRUCTURE_H
