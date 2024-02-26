// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME BIBDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "BIB.hh"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *BIB_Dictionary();
   static void BIB_TClassManip(TClass*);
   static void *new_BIB(void *p = nullptr);
   static void *newArray_BIB(Long_t size, void *p);
   static void delete_BIB(void *p);
   static void deleteArray_BIB(void *p);
   static void destruct_BIB(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::BIB*)
   {
      ::BIB *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::BIB));
      static ::ROOT::TGenericClassInfo 
         instance("BIB", "BIB.hh", 17,
                  typeid(::BIB), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &BIB_Dictionary, isa_proxy, 4,
                  sizeof(::BIB) );
      instance.SetNew(&new_BIB);
      instance.SetNewArray(&newArray_BIB);
      instance.SetDelete(&delete_BIB);
      instance.SetDeleteArray(&deleteArray_BIB);
      instance.SetDestructor(&destruct_BIB);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::BIB*)
   {
      return GenerateInitInstanceLocal(static_cast<::BIB*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::BIB*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *BIB_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::BIB*>(nullptr))->GetClass();
      BIB_TClassManip(theClass);
   return theClass;
   }

   static void BIB_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_BIB(void *p) {
      return  p ? new(p) ::BIB : new ::BIB;
   }
   static void *newArray_BIB(Long_t nElements, void *p) {
      return p ? new(p) ::BIB[nElements] : new ::BIB[nElements];
   }
   // Wrapper around operator delete
   static void delete_BIB(void *p) {
      delete (static_cast<::BIB*>(p));
   }
   static void deleteArray_BIB(void *p) {
      delete [] (static_cast<::BIB*>(p));
   }
   static void destruct_BIB(void *p) {
      typedef ::BIB current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::BIB

namespace {
  void TriggerDictionaryInitialization_BIBDict_Impl() {
    static const char* headers[] = {
"BIB.hh",
nullptr
    };
    static const char* includePaths[] = {
"/lustre/cmswork/fnardi/anaconda3/envs/geant4/include/",
"/home/fnardi/bib_tuples/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "BIBDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$BIB.hh")))  BIB;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "BIBDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "BIB.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"BIB", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("BIBDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_BIBDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_BIBDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_BIBDict() {
  TriggerDictionaryInitialization_BIBDict_Impl();
}
