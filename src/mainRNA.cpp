
#include"rnastructure.h"
#include<iostream>
#include"fstream.h"


using namespace std;



int main(){
    ifstream input;
    input.open("/Users/margrethmpossi/Desktop/Previous classes/GENE 218/RNAStruct/res/tRNA.txt");
    RNAStructure stru(input);
    cout << "Please enter a starting value:" << endl;
    return 0;
}
