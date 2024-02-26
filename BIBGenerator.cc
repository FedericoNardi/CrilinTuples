// BIBGenerator.cc
#include "BIB.hh"

int main(){
   BIB* bib;
   for(int ii=0; ii<12; ii++){
      bib->Loop(ii, false, false, true);
   }
   delete bib;
   return 0;
}