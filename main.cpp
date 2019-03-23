
#include <stdlib.h>
#include "ProtocolParty.h"
#include "BAParty.h"
#include "ZpKaratsubaElement.h"
#include <smmintrin.h>
#include <inttypes.h>
#include <stdio.h>
#include <x86intrin.h>


template <class FieldType>
void testReconstruct(const vector<int>& poly,
                     const vector<int>& alpha,
                     vector<int>& re_poly){
  int nCoeff = poly.size();
  int nPoints = alpha.size();
  vector<FieldType> fieldPoly(nCoeff, FieldType());
  vector<FieldType> fieldX(nPoints, FieldType());
  vector<FieldType> fieldY(nPoints, FieldType());
  vector<FieldType> result;

  // build polynomial
  for(int i=0; i<nCoeff; i++){
    fieldPoly[i]= FieldType( poly[i] );
  }

  // build sample points
  for(int i=0; i<nPoints; i++){
    fieldX[i] = FieldType( alpha[i] );
    fieldY[i] = evalPolynomial<FieldType>( fieldPoly, fieldX[i] );
    
  }

  // try reconstruction
  // reconstruct <FieldType> (fieldX, fieldY, nCoeff-1, result);

  cout << "evaluation result:" << endl;
  
  for(int i=0; i<nPoints; i++){
    cout << "(" << fieldX[i] << " " << fieldY[i] << ")" << endl;
  }
  return;
}


int main(int argc, char* argv[])
{

    CmdParser parser;
    auto parameters = parser.parseArguments("", argc, argv);
    int times = stoi(parser.getValueByKey(parameters, "internalIterationsNumber"));


    string fieldType = parser.getValueByKey(parameters, "fieldType");
    cout<<"fieldType = "<<fieldType<<endl;

    int polynomial[] = {4, 1, 2, 3, 5}; // size = 5
    int alpha[] = {1, 4, 9, 20, 5, 8, 6, 3, 7, 2}; // size = 10
    
    vector<int> result;
    vector<int> poly(polynomial, polynomial+5);
    vector<int> alph(alpha, alpha+10);

    if(fieldType.compare("ZpMersenne31") == 0)
      {
        testReconstruct<ZpMersenneIntElement>(poly, alph, result);
      }
    else if(fieldType.compare("ZpMersenne61") == 0)
      {
        testReconstruct<ZpMersenneLongElement>(poly, alph, result);
      }
    else if(fieldType.compare("ZpKaratsuba") == 0)
      {
        testReconstruct<ZpKaratsubaElement>(poly, alph, result);
      }
    else if(fieldType.compare("GF2m") == 0)
      {
        testReconstruct<GF2E>(poly, alph, result);
      }
    else if(fieldType.compare("Zp") == 0)
      {
        testReconstruct<ZZ_p>(poly, alph, result);
      }

    return 0;
}
