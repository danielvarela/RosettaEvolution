
#ifndef types_hpp
#define types_hpp
#include <string>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <iomanip> /* std::setprecision */
#include <algorithm>    // std::sort
#include <fstream>
#include <random>
#include <cmath>

using namespace std;

const double PI = std::atan(1.0)*4;
const double RAD2DEG = 180/PI;
const double DEG2RAD = PI/180;
const int NUM_DAT = 6;
const int NUM_MIS = 5;

// amino acid residue types (20 classes)
const int ALA = 0;
const int CYS = 1;
const int ASP = 2;
const int GLU = 3;
const int PHE = 4;
const int GLY = 5;
const int HIS = 6;
const int ILE = 7;
const int LYS = 8;
const int LEU = 9;
const int MET = 10;
const int ASN = 11;
const int PRO = 12;
const int GLN = 13;
const int ARG = 14;
const int SER = 15;
const int THR = 16;
const int VAL = 17;
const int TRP = 18;
const int TYR = 19;


// secondary structure types
const int HELIX = 0;
const int SHEET = 1;
const int COIL = 2;
const int UNK = 3;


// structure for storing fragment information
struct FRAGMENT {
  int POS;
  int FRG;
  int ID;
  string AA;
  string SS;
  double PHI;
  double PSI;
  double OME;
};

//bool compareEntry( const FRAGMENT &lhs, const FRAGMENT &rhs);
/*************************************************************************
 * Name        : compareEntry
 * Purpose     : helper function for multi-column sorting
 * Arguments   : const FRAGMENT &lhs, const FRAGMENT &rhs
 * Return Type : bool
 *************************************************************************/
inline bool compareEntry( const FRAGMENT &lhs, const FRAGMENT &rhs ) {
  if (lhs.POS < rhs.POS)
    return true;
  else if (rhs.POS < lhs.POS)
    return false;
  else // POS equals.
    if (lhs. FRG  < rhs. FRG )
      return true;
    else if (rhs. FRG  < lhs. FRG )
      return false;
    else // POS and FRG equals
      return lhs.ID < rhs.ID;
}


/*************************************************************************
 * Name        : writeFragment
 * Purpose     : writes fragment in ROSETTA's format
 * Arguments   : vector<FRAGMENT> & fragment,ostream & out
 * Return Type : void
 *************************************************************************/
inline void writeFragment(vector<FRAGMENT> & fragment,ostream & out, int neighbors_ = 2) {
    int prev = -1;
    int prev_frg = -1;
    int F = fragment.size();
    for (int i =0 ; i < F; i++) {
        if (fragment[i].POS != prev && prev == -1) {
            out <<  " position:          "<<setw(3)<< std::right<<fragment[i].POS<<" neighbors:          "<<setw(3)<< std::right<<neighbors_<<endl;
            prev = fragment[i].POS;
            if (fragment[i].FRG != prev_frg) {
                prev_frg = fragment[i].FRG;
                out << endl;
                out << " xxxx A" << setw(6)<< std::right<<fragment[i].ID << " " << fragment[i].AA << " " <<fragment[i].SS <<
                " " <<setw(8) <<fixed << showpoint <<std::right <<std::setprecision(3)<<fragment[i].PHI<<
                " " <<setw(8) <<fixed << showpoint <<std::right <<std::setprecision(3)<<fragment[i].PSI<<
                " " <<setw(8) <<fixed << showpoint <<std::right <<std::setprecision(3)<<fragment[i].OME<<
                "    0.000    0.000     0.000 1     0.000 P"<<setw(3)<< std::right<<fragment[i].POS<<" F"<<
                setw(3)<< std::right<<fragment[i].FRG <<endl;
            }
            else {
                out << endl;
                out << " xxxx A" << setw(6)<< std::right<<fragment[i].ID << " " << fragment[i].AA << " " <<fragment[i].SS <<
                " " <<setw(8) <<fixed << showpoint <<std::right <<std::setprecision(3)<<fragment[i].PHI<<
                " " <<setw(8) <<fixed << showpoint <<std::right <<std::setprecision(3)<<fragment[i].PSI<<
                " " <<setw(8) <<fixed << showpoint <<std::right <<std::setprecision(3)<<fragment[i].OME<<
                "    0.000    0.000     0.000 1     0.000 P"<<setw(3)<< std::right<<fragment[i].POS<<" F"<<
                setw(3)<< std::right<<fragment[i].FRG <<endl;
                
            }
        }
        else if (fragment[i].POS != prev && prev != -1) {
            out << endl;
            out <<  " position:          "<<setw(3)<< std::right<<fragment[i].POS<<" neighbors:          "<<setw(3)<< std::right<<neighbors_<<endl;
            prev = fragment[i].POS;
            if (fragment[i].FRG != prev_frg) {
                prev_frg = fragment[i].FRG;
                out << endl;
                out << " xxxx A" << setw(6)<< std::right<<fragment[i].ID << " " << fragment[i].AA << " " <<fragment[i].SS <<
                " " <<setw(8) <<fixed << showpoint <<std::right <<std::setprecision(3)<<fragment[i].PHI<<
                " " <<setw(8) <<fixed << showpoint <<std::right <<std::setprecision(3)<<fragment[i].PSI<<
                " " <<setw(8) <<fixed << showpoint <<std::right <<std::setprecision(3)<<fragment[i].OME<<
                "    0.000    0.000     0.000 1     0.000 P"<<setw(3)<< std::right<<fragment[i].POS<<" F"<<
                setw(3)<< std::right<<fragment[i].FRG <<endl;
            }
            else {
                out << endl;
                out << " xxxx A" << setw(6)<< std::right<<fragment[i].ID << " " << fragment[i].AA << " " <<fragment[i].SS <<
                " " <<setw(8) <<fixed << showpoint <<std::right <<std::setprecision(3)<<fragment[i].PHI<<
                " " <<setw(8) <<fixed << showpoint <<std::right <<std::setprecision(3)<<fragment[i].PSI<<
                " " <<setw(8) <<fixed << showpoint <<std::right <<std::setprecision(3)<<fragment[i].OME<<
                "    0.000    0.000     0.000 1     0.000 P"<<setw(3)<< std::right<<fragment[i].POS<<" F"<<
                setw(3)<< std::right<<fragment[i].FRG <<endl;
                
            }
        }
        else {
            if (fragment[i].FRG != prev_frg) {
                prev_frg = fragment[i].FRG;
                out << endl;
                out << " xxxx A" << setw(6)<< std::right<<fragment[i].ID << " " << fragment[i].AA << " " <<fragment[i].SS <<
                " " <<setw(8) <<fixed << showpoint <<std::right <<std::setprecision(3)<<fragment[i].PHI<<
                " " <<setw(8) <<fixed << showpoint <<std::right <<std::setprecision(3)<<fragment[i].PSI<<
                " " <<setw(8) <<fixed << showpoint <<std::right <<std::setprecision(3)<<fragment[i].OME<<
                "    0.000    0.000     0.000 1     0.000 P"<<setw(3)<< std::right<<fragment[i].POS<<" F"<<
                setw(3)<< std::right<<fragment[i].FRG <<endl;
            }
            else {
                
                out << " xxxx A" << setw(6)<< std::right<<fragment[i].ID << " " << fragment[i].AA << " " <<fragment[i].SS <<
                " " <<setw(8) <<fixed << showpoint <<std::right <<std::setprecision(3)<<fragment[i].PHI<<
                " " <<setw(8) <<fixed << showpoint <<std::right <<std::setprecision(3)<<fragment[i].PSI<<
                " " <<setw(8) <<fixed << showpoint <<std::right <<std::setprecision(3)<<fragment[i].OME<<
                "    0.000    0.000     0.000 1     0.000 P"<<setw(3)<< std::right<<fragment[i].POS<<" F"<<
                setw(3)<< std::right<<fragment[i].FRG <<endl;

            }
        }
    }
}


#endif // types_hpp
