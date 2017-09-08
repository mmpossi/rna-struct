#include "rnastructure.h"
#include<string>
#include "Vector.h"
#include"Set.h"
#include <iostream>
#include <istream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <sstream>
#include "SimpleGraph.h"
#include <map.h>
#include <queue>
#include <iterator>


// ALTERNATIVE: use the in-loop probabilities derived from count data to calculate total probability
// total probability is equal to  sum of log likelihoods for 3-7 neighbouring nucleotides.

//Defining max
//storing structure only if larger than max  set<pair<int, int> > max;
//update max as you go   if (copyPair0> basePairs) ->if (copyPair0> max){
//saving only the modified structures,not the original structure. Adding empty after each iteration
//removing wobble bps
//iterating through available in symmetry and countBase
//num =20 ARBITRARILY chosen as the number limit for storing duplicates and starting empyty makes
//adding edit model save scaffold model NOT YET
//Changed evaluation of countBase return, should strictly be larger than 2. Rambled structures removed
//change structure of m_interests to vector. Set is strictly ordereed.




using namespace std;

bool RNAStructure::compare(set<pair<int, int> > a, set<pair<int, int> > b){
    return a.size() < b.size();
}


RNAStructure::RNAStructure(ifstream& input){
    generate(input);
}


void RNAStructure:: generate(ifstream& input){
    vector<int>  intePositions;// position with sequences likely to be loops
    initiate(input,intePositions);
    m_interests=intePositions;
    set<int> available =m_available;// all unpaired position
    set<pair<int, int> >  empty;
    //int order=0;
    int prevPosition=0;
    queue <set<pair<int, int> >> structures;// all possible structures
    structures.push(empty);// includes empty in the list of structures
    set<pair<int, int> > max;//Defining max
    vector<int> ::iterator it1=(m_interests).begin();
    while (it1!=m_interests.end()) {

        int num=structures.size();
        int index=*it1;



        //set<pair<int, int> >  basePairs;
        if (num<m_length/10) {
            set<pair<int, int> > copyPair0;
            copyPair0=coinStructure(available, empty, index);
            structures.push(copyPair0);}
        for (size_t i=0; i<num; i++ ){//loops through all previous structures and modifies using new intePosition

            set<pair<int, int> >  basePairs;
            set<pair<int, int> > copyPair0;
            basePairs=structures.front();
            //structures.pushback(basePairs);
            structures.pop();
            set<int> available =m_available;
            updateAvailable(basePairs,available);
            copyPair0=coinStructure(available, basePairs, index);

            cout<<"Num of base pairs is "<< copyPair0.size()<<endl;
            if ( (copyPair0.size()-basePairs.size())>2) {
                structures.push(copyPair0);//checks if new structure is more favourable that starting structure
                if (copyPair0.size() > max.size()){
                    max=copyPair0;
                    cout<<"Num of shown base pairs is "<< copyPair0.size()<<endl;
                    updateStruct(copyPair0);
                }
            }
            else {
                structures.push(basePairs);
            }
        }
        structures.push(empty);
        prevPosition=index;
        it1++;
    }
    cout << "Total num of Structures is " <<structures.size()<< endl;
    /* do {

        empty=structures.front();
        cout<<"Num of base pairs is "<< empty.size()<<endl;
        updateStruct(empty);
        structures.pop();
    } while (!structures.empty() );*/


}


void RNAStructure::initiate (ifstream& input, vector<int> & intePositions){

    for (int i=0; i<m_loopPair.size(); i++){
        vector<int> values;
        intePositions=(values);//inserts list of interesting positions by the dinucleotide frequency in loops
    }

    char A='A';
    char T='T';
    char U='U';
    char G='G';
    char C='C';
    vector<char> Acomplements;
    Acomplements.push_back(U);
    Acomplements.push_back(T);
    m_antisense[A]= Acomplements;

    vector<char> Gcomplements;
    Gcomplements.push_back(C);
    Gcomplements.push_back(U);
    Gcomplements.push_back(T);
    m_antisense[G]= Gcomplements;

    vector<char> Tcomplements;
    Tcomplements.push_back(A);
    Tcomplements.push_back(G);
    m_antisense[T]=Tcomplements;

    vector<char> Ucomplements;
    Ucomplements.push_back(A);
    Ucomplements.push_back(G);
    m_antisense[U]= Ucomplements;

    vector<char> Ccomplements;
    Ccomplements.push_back(G);
    m_antisense[C] =Ccomplements;
    Ccomplements.clear();
    Tcomplements.clear();
    Ucomplements.clear();
    Acomplements.clear();
    Gcomplements.clear();
    string line;
    string twos;
    while (!input.eof()){
        getline(input, line);
        //cout << line<< endl;
        int len=line.length();
        for (int i=0; i< len; i++){
            m_sequence.push_back(line[i]);
            cout << line[i] << "is read"<< endl;
            int num=m_sequence.size();
            m_available.insert(num-1);
            twos=line.substr(i,2);
        }
    }

    /**Reading the file from the classifier
     */
    ifstream inputPositions;
    inputPositions.open("/Users/margrethmpossi/Desktop/Previous classes/GENE 218/RNAStruct/res/positions");
    string line2;
    getline(inputPositions, line2);
    while(line2!=""){
        istringstream lin(line2);
        cout << "line read is "<<line2 << endl;
        int token;
        lin >>token;
        while (lin>>token){
            intePositions.push_back(token);
        }
        getline(inputPositions, line2);
    }
    m_numInterests=intePositions.size();//stores the number of positions of interest
    m_length=m_sequence.size()-1;
    initiateGraph();
}


/**
 * @brief RNAStructure::initiateGraph
 *copies the m_sequence into a graph
 */
void RNAStructure::initiateGraph(){
    inCircle(m_length+3,graph);
}

//sets up the nodes vector in graph
void RNAStructure:: inCircle(const int num, SimpleGraph& graph){
    const double pi= 3.14159265358979323;
    int i=0;
    vector<char>::iterator it=m_sequence.begin();
    while (it!=m_sequence.end()){
        Node node;
        Edge edge;
        double xValue = (cos((2*pi*i)/num));
        double yValue = (sin((2*pi*i)/num));
        node.x=xValue;
        node.y=yValue;
        node.type=*it;
        graph.nodes.push_back(node);
        if (i==0 ){
            edge.start=i;}
        else {
            if (i==num){
                edge.end=edge.start;
                edge.start=0;
                graph.edges.push_back(edge);
            }else{
                edge.end=edge.start;
                edge.start=i;
                graph.edges.push_back(edge);
            }
        }
        it++;
        i++;
    }

    DrawGraph(graph);
}

set<pair<int, int> >  RNAStructure:: coinStructure( set<int> available,  set<pair<int, int> > basePairs ,  int index){
    symmetry(available, index, basePairs); //should return optimal solution
    return basePairs;
}

void RNAStructure:: updateAvailable(set<pair<int, int> >& basePairs,set<int>& available ){
    set <pair<int, int> >::iterator it=basePairs.begin();
    while (it!=basePairs.end()){
        pair <int, int> bond=*it;
        set<int>::iterator it2=available.find(min(bond.first, bond.second));
        if (it2!=available.end()){
            available.erase(it2);
            it2=available.find(max(bond.first, bond.second));
            if (it2!=available.end()){
                available.erase(it2);
            }
        }
        it++;
    }
}



void RNAStructure:: symmetry( set<int> available,int index,  set<pair<int, int> >& basePairs){//EDITED
    set <pair<int, int> > copyPairs=basePairs;//previous pairs
    set <pair<int, int> > prevPairs=basePairs;
    //  if (index==0||index==m_length) return;

    set<int>::iterator it1=available.find(index);//checks downstream
    set<int>::iterator it2=available.find(index);//checks upstream

  /*  for (int i=0; i<RANGE; i++){
        if (it1!=available.begin())it1--;
        if (it2!=available.end())it2++;
    }*/
     for (int i=0; i<RANGE; i++){
        vector <char> complements=m_antisense.at(m_sequence[*it1]);
         for (int j=0; i<RANGE; j++){
            // if ((it2>it1&&(it2-it1)>(index-*it1) )||( it2<it1&& (it2-it1)<(index-*it1))){
            if (toCheck(*it2, complements, available)){
                copyPairs=countBase(it1,it2, basePairs, available);
                if(copyPairs.size()>prevPairs.size()+2)  {
                    prevPairs=copyPairs;
                }
            }
            // }
            if (it2!=available.end()) it2++; else {break;}
        }
        if (it1!=available.begin())it1--; else {break;}
    }
    if(prevPairs.size()>basePairs.size())  {
        basePairs=prevPairs;
    }
    //while ...




}


/*void RNAStructure:: symmetry( set<int> available,int index,  set<pair<int, int> >& basePairs){//TO DO
    set <pair<int, int> > copyPairs=basePairs;//previous pairs
    set <pair<int, int> > prevPairs=basePairs;
    if (index==0||index==m_length) return;
    else{
        for(int ai=max(index-RANGE,0); ai<min(index+RANGE, m_length); ai++){
            set<int>::iterator it=available.find(ai);
            if(it!=available.end()){
                vector <char> complements=m_antisense.at(m_sequence[ai]);
                for (int bij= max(ai-RANGE,0); bij<min((ai+RANGE), m_length); bij++){
                    if ((bij>ai&&(bij-ai)>(index-ai) )||( bij<ai&& (bij-ai)<(index-ai))){
                        if (toCheck(bij, complements, available)){
                            copyPairs=countBase(min(ai,bij), max(ai, bij), basePairs, available);
                            if(copyPairs.size()>prevPairs.size())  {
                                prevPairs=copyPairs;
                            }
                        }
                    }
                }}
        }
        if(prevPairs.size()>basePairs.size())  {
            basePairs=prevPairs;
        }
    }
}*/


set<pair<int, int> >  RNAStructure::countBase(set<int>::iterator it1,set<int>::iterator it2,  set <pair <int, int> > basePairs, set<int>  available){
    pair<int, int> bond;
    bond.first=*it1;
    bond.second=*it2;
    //cout<<"changing iterator"<< endl;
    basePairs.insert(bond);
    if (it1!=available.end() && it2!=available.begin()){
        available.erase(*it1);
        available.erase(*it2);
        --it2;
        ++it1;

    }
    cout<<*it1<< " , "<<*it2 <<endl;

    if (it1!=available.end() && it2!=available.begin()) {
        //set<int>::iterator it=available.find(start-1);
        //if (it!=available.end()){
        vector<char> complements=m_antisense[m_sequence[*it1]];
        if (toCheck (*it2, complements, available)) {
            basePairs=countBase(it1, it2, basePairs,available);
            return basePairs;
        }
    }
    return basePairs;
}

/*
set<pair<int, int> >  RNAStructure::countBase(int start,int end,  set <pair <int, int> > basePairs, set<int>  available){
    pair<int, int> bond;
    bond.first=start;
    bond.second=end;
    basePairs.insert(bond);
    available.erase(available.find(start));
    available.erase(available.find(end));
    if (start-1>=0 && end+1<=m_length) {
        set<int>::iterator it=available.find(start-1);
        if (it!=available.end()){
            vector<char> complements=m_antisense[m_sequence[start-1]];
            if (toCheck (end+1, complements, available)) {
                basePairs=countBase(start-1, end+1, basePairs,available);
                return basePairs;
            }}
    }
    return basePairs;
}
*/

/**
 * @brief RNAStructure::toCheck
 * @param index
 * @param complements
 * @return
 *checks if the index should be checked for potential pairing with another
 */
bool RNAStructure::toCheck(int index, vector <char>& complements,const set<int> & available){
    if (index<0|| index>m_length) return false;
    char ch=m_sequence[index];
    vector<char>::iterator it=complements.begin();
    if (it!=complements.end()&& *it==ch){
        set<int>::iterator it2=available.find(index);
        if (it2!=available.end()){
            return true;
        }
    }
    return false;
}



void RNAStructure:: updateStruct(const set <pair<int, int> >& basePairs){
    SimpleGraph graphCopy=graph;
    cout<<"Num of shown base pairs is "<< basePairs.size()<<endl;
    set <pair<int, int> >::iterator it=basePairs.begin();
    while (it!=basePairs.end()){
        graphCopy.nodes.at((*it).first).bond=(*it).second;
        graphCopy.nodes.at((*it).second).bond=(*it).first;
        Edge edge;
        edge.start=(*it).first;
        edge.end=(*it).second;
        graphCopy.edges.push_back(edge);\
        it++;
    }
    sleep(1);
    DrawGraph(graphCopy);
}


//return a vector for relative to node1
pair <double, double> RNAStructure:: getVector(Node node1, Node node2){
    pair<double, double> position;
    double x1=node1.x;
    double x2=node2.x;
    double y1=node1.y;
    double y2=node2.y;
    position.first=sqrt (pow((x1 -x2), 2.0) +pow(( y1 - y2), 2.0));
    position.second=atan2(y1-y2,x1-x2);
    return position;
}
// takes a pair of nodes and returns the force between them
pair<double, double>RNAStructure:: getForce(Node node1, Node node2,double KAttract){
    pair<double,double> force;
    pair<double,double> position=getVector(node1,node2);
    double forceTotal=(kRepel/position.first) -(KAttract/position.first);
    force.first =forceTotal*cos(position.second);
    force.second =forceTotal*sin(position.second);

    return force;

}
//returns total distance change after all the force
pair<double, double>RNAStructure:: totalForce(Node& node1, Node& node2, Node& node3, Node& node4){
    pair<double, double> totalForce;
    pair<double, double> totalForce1;
    pair<double, double> totalForce2;
    totalForce= getForce(node1, node2,kAttract);
    totalForce1=getForce(node1, node3,kAttract);
    totalForce2=getForce(node1,node4, 10.0*kAttract);
    double delX=totalForce.first;
    double delY=totalForce.second;
    node2.x=node2.x-delX;
    node2.y=node2.y-delY;

    delX=totalForce1.first;
    delY=totalForce1.second;
    node3.x=node3.x-delX;
    node3.y=node3.y-delY;

    delX=totalForce2.first;
    delY=totalForce2.second;
    node4.x=node4.x-delX;
    node4.y=node4.y-delY;

    double x =totalForce.first+totalForce1.first+totalForce2.first;
    double y =totalForce.second+totalForce1.second+totalForce2.second;

    node1.x=node1.y+ x;
    node1.y=node1.y+y ;
    return totalForce;
}


void RNAStructure:: optimize(SimpleGraph& graphy, const set<pair<int, int> >& basePairs){
    pair<double, double> force;
    set <pair<int, int> >::iterator it=basePairs.begin();
    while (it!=basePairs.end()){
        Node node1=graphy.nodes.at((*it).first);
        Node node2;
        Node node3;
        Node node4=graphy.nodes.at((*it).second);
        if(!(*it).first==0){
            // cout >> (*it).first >> " na " >> (*it).second>> endl;//testing base pair arrangement
            node2=graphy.nodes.at((*it).first-1);
            if((*it).first!=m_length){
                node3=graphy.nodes.at((*it).first+1);
                force =totalForce(node1, node2, node3, node4);
            }
        }
        it++;
    }
}

