#ifndef RNASTRUCTURE_H
#define RNASTRUCTURE_H

/*This class was written by Margreth Mpossi.
 *The program used simple graph to represent RNA secondary structure.
 *Functions in this class work together to produce low-energy stem-loops
  */


#include<iostream>
#include <fstream.h>
#include <sstream>
#include<string>
#include "Vector.h"
#include"Set.h"
#include"map.h"
#include "SimpleGraph.h"
#include <queue>



class RNAStructure
{
public:
    /*
     *Initiates RNA structure from an ifstream reader
     */
    RNAStructure(std::ifstream& input);



private:
    /*
     *The range of indeces to search when symmetry() is called
     */
    const int RANGE=5;

    /*
     *The minimum number of base pairs considered significant
     * @brief MIN_PAIRS
     */
    const int MIN_PAIRS=3;

    /*
     * @brief m_numInterests
     *Total number of all indices in the m_interest vector.
     *The number of indices to call symmetry for
     */
    int m_numInterests;

    /*
     *Number of pairs considered significant determiners of loop locations
     * @brief m_loopPair
     */
    vector<std::string> m_loopPair;// all sequences frequently found in loops

    /*
     *A list of all indices with dinucleotides of interest
     * @brief m_interests
     */
    vector<int> m_interests;// stores list of indeces of interest

    /*List of all initially non-paired positions
     *Essentially the entire sequence
     * @brief m_available
     */
    set <int>m_available;

    /*the entire sequence as a vector, for direct access of an index
     * @brief m_sequence
     */
    vector< char>  m_sequence;//all nucleotides in the sequence, position same as the index in the vector

    /* a map of all bases and their complementary W/C pairs
     * @brief m_antisense
     */
    map<char, vector<char >> m_antisense;// stores values and their antisense for easy location of complements

    std::queue <set<pair<int, int> >> m_structures;

    /* length of the sequence
     * @brief m_length
     */
    int m_length;//length of the sequence

    /*a graph of all unpairs nucleotides in a circle
     * @brief graph
     */

    SimpleGraph graph;
    /*
     *the time to keep iterating for distance optimization
     * @brief m_time
     */
    const double m_time=70;

    /* contant for repulsion when optiizing the secondary sctructure
     * @brief kRepel
     */
    const double kRepel=0.0009;

    /*constant for attraction when optinizing represantation of final structure
     * @brief kAttract
     */
    const double kAttract=0.001;

    /*Inititiates RNA structure, including the sequence and the positions of interest
     *Positions of interest are noted while reading the input file
     * @brief initiate
     * @param input
     * @param intePositions
     */
    void initiate (ifstream& input, vector<int >& intePositions);

    /*iterates through positions of interest and returns the values of
     * @brief generate
     * @param input
     */
    void generate(std::ifstream& input);

    /*Initiates the graph in a circle
     * @brief initiateGraph
     */
    void initiateGraph();// copies m_sequence into a graph

    /*represents each position as a node in a graph in a circular conformation
     * @brief inCircle
     * @param num
     * @param graph
     */
    void inCircle(const int num, SimpleGraph& graph); //initializes positions in a circle

    /*Looks for complementarity around the given indices from among the available positions
     *adds findings to the list of base pairs
     * @brief symmetry
     * @param available
     * @param index
     * @param basePairs
     */
    void  symmetry( set<int> available,int index,  set<pair<int, int> >& basePairs);// for each position of interest, finds symmetry

    /*returns if an index is available and a acomplement of a given partner position
     * @brief toCheck
     * @param index
     * @param complements
     * @param available
     * @return
     */
    bool toCheck(int index, vector <char>& complements,const set<int> & available);

    /*
     * @brief countBase
     * @param start
     * @param end
     * @param basePairs
     * @param available
     * @return
     */
    set <pair <int, int> > countBase(set<int>::iterator it1,set<int>::iterator it2,  set <pair <int, int> > basePairs, set<int>  available);

    /* updates graph display based on the base pairs present
     * @brief updateStruct
     * @param basePairs
     * @param available
     */
    void updateStruct(const set <pair<int, int> >& basePairs);

    /*returns base pairs from optimizing symmetry at an index
      * @brief coinStructure
      * @param available
      * @param basePairs
      * @param index
      * @return
      */
     set<pair<int, int> >  coinStructure( set<int> available,  set<pair<int, int> > basePairs ,  int index);
   /*comparator function for the priority queue
    */
     bool compare(set<pair<int, int> > a, set<pair<int, int> > b);

     /*updates the list of available positions, by removing all positions in the set of base pairs
     * @brief updateAvailable
     * @param basePairs
     * @param available
     */
    void updateAvailable(set<pair<int, int> >& basePairs,set<int>& available );

    /*calculates the direction vector between to nodes
     * @brief getVector
     * @param node1
     * @param node2
     * @return
     */
    pair <double, double> getVector(Node node1, Node node2);
    /*Calculats the force betwee two nodes from the direction vector and the force constants, kRepel and kAttract
     * @brief getForce
     * @param node1
     * @param node2
     * @param KAttract
     * @return
     */
    pair<double, double> getForce(Node node1, Node node2,double KAttract);

    /*Calculates the change in distance between two nodes from force calculations
     * @brief totalForce
     * @param node1
     * @param node2
     * @param node3
     * @param node4
     * @return
     */
     pair<double, double> totalForce(Node& node1, Node& node2, Node& node3, Node& node4);
    /* updates position of nodes from distance changes in force calculations
     * @brief optimize
     * @param graphy
     * @param basePairs
     */
    void optimize(SimpleGraph& graphy, const set<pair<int, int> >& basePairs);

    void reposition(SimpleGraph& graphy);



};



#endif // RNASTRUCTURE_H
