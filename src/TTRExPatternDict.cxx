// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TTRExPatternDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
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

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "TTRExPattern.hxx"
#include "OutLink.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *trexcLcLTTRExEvent_Dictionary();
   static void trexcLcLTTRExEvent_TClassManip(TClass*);
   static void *new_trexcLcLTTRExEvent(void *p = 0);
   static void *newArray_trexcLcLTTRExEvent(Long_t size, void *p);
   static void delete_trexcLcLTTRExEvent(void *p);
   static void deleteArray_trexcLcLTTRExEvent(void *p);
   static void destruct_trexcLcLTTRExEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::trex::TTRExEvent*)
   {
      ::trex::TTRExEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::trex::TTRExEvent));
      static ::ROOT::TGenericClassInfo 
         instance("trex::TTRExEvent", "TTRExPattern.hxx", 14,
                  typeid(::trex::TTRExEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &trexcLcLTTRExEvent_Dictionary, isa_proxy, 4,
                  sizeof(::trex::TTRExEvent) );
      instance.SetNew(&new_trexcLcLTTRExEvent);
      instance.SetNewArray(&newArray_trexcLcLTTRExEvent);
      instance.SetDelete(&delete_trexcLcLTTRExEvent);
      instance.SetDeleteArray(&deleteArray_trexcLcLTTRExEvent);
      instance.SetDestructor(&destruct_trexcLcLTTRExEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::trex::TTRExEvent*)
   {
      return GenerateInitInstanceLocal((::trex::TTRExEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::trex::TTRExEvent*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *trexcLcLTTRExEvent_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::trex::TTRExEvent*)0x0)->GetClass();
      trexcLcLTTRExEvent_TClassManip(theClass);
   return theClass;
   }

   static void trexcLcLTTRExEvent_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *trexcLcLTTRExPattern_Dictionary();
   static void trexcLcLTTRExPattern_TClassManip(TClass*);
   static void *new_trexcLcLTTRExPattern(void *p = 0);
   static void *newArray_trexcLcLTTRExPattern(Long_t size, void *p);
   static void delete_trexcLcLTTRExPattern(void *p);
   static void deleteArray_trexcLcLTTRExPattern(void *p);
   static void destruct_trexcLcLTTRExPattern(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::trex::TTRExPattern*)
   {
      ::trex::TTRExPattern *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::trex::TTRExPattern));
      static ::ROOT::TGenericClassInfo 
         instance("trex::TTRExPattern", "TTRExPattern.hxx", 49,
                  typeid(::trex::TTRExPattern), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &trexcLcLTTRExPattern_Dictionary, isa_proxy, 4,
                  sizeof(::trex::TTRExPattern) );
      instance.SetNew(&new_trexcLcLTTRExPattern);
      instance.SetNewArray(&newArray_trexcLcLTTRExPattern);
      instance.SetDelete(&delete_trexcLcLTTRExPattern);
      instance.SetDeleteArray(&deleteArray_trexcLcLTTRExPattern);
      instance.SetDestructor(&destruct_trexcLcLTTRExPattern);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::trex::TTRExPattern*)
   {
      return GenerateInitInstanceLocal((::trex::TTRExPattern*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::trex::TTRExPattern*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *trexcLcLTTRExPattern_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::trex::TTRExPattern*)0x0)->GetClass();
      trexcLcLTTRExPattern_TClassManip(theClass);
   return theClass;
   }

   static void trexcLcLTTRExPattern_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *trexcLcLTTPCHitPad_Dictionary();
   static void trexcLcLTTPCHitPad_TClassManip(TClass*);
   static void *new_trexcLcLTTPCHitPad(void *p = 0);
   static void *newArray_trexcLcLTTPCHitPad(Long_t size, void *p);
   static void delete_trexcLcLTTPCHitPad(void *p);
   static void deleteArray_trexcLcLTTPCHitPad(void *p);
   static void destruct_trexcLcLTTPCHitPad(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::trex::TTPCHitPad*)
   {
      ::trex::TTPCHitPad *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::trex::TTPCHitPad));
      static ::ROOT::TGenericClassInfo 
         instance("trex::TTPCHitPad", "TTPCHitPad.hxx", 18,
                  typeid(::trex::TTPCHitPad), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &trexcLcLTTPCHitPad_Dictionary, isa_proxy, 4,
                  sizeof(::trex::TTPCHitPad) );
      instance.SetNew(&new_trexcLcLTTPCHitPad);
      instance.SetNewArray(&newArray_trexcLcLTTPCHitPad);
      instance.SetDelete(&delete_trexcLcLTTPCHitPad);
      instance.SetDeleteArray(&deleteArray_trexcLcLTTPCHitPad);
      instance.SetDestructor(&destruct_trexcLcLTTPCHitPad);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::trex::TTPCHitPad*)
   {
      return GenerateInitInstanceLocal((::trex::TTPCHitPad*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::trex::TTPCHitPad*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *trexcLcLTTPCHitPad_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::trex::TTPCHitPad*)0x0)->GetClass();
      trexcLcLTTPCHitPad_TClassManip(theClass);
   return theClass;
   }

   static void trexcLcLTTPCHitPad_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *trexcLcLTTRExHVCluster_Dictionary();
   static void trexcLcLTTRExHVCluster_TClassManip(TClass*);
   static void *new_trexcLcLTTRExHVCluster(void *p = 0);
   static void *newArray_trexcLcLTTRExHVCluster(Long_t size, void *p);
   static void delete_trexcLcLTTRExHVCluster(void *p);
   static void deleteArray_trexcLcLTTRExHVCluster(void *p);
   static void destruct_trexcLcLTTRExHVCluster(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::trex::TTRExHVCluster*)
   {
      ::trex::TTRExHVCluster *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::trex::TTRExHVCluster));
      static ::ROOT::TGenericClassInfo 
         instance("trex::TTRExHVCluster", "TTRExHVCluster.hxx", 16,
                  typeid(::trex::TTRExHVCluster), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &trexcLcLTTRExHVCluster_Dictionary, isa_proxy, 4,
                  sizeof(::trex::TTRExHVCluster) );
      instance.SetNew(&new_trexcLcLTTRExHVCluster);
      instance.SetNewArray(&newArray_trexcLcLTTRExHVCluster);
      instance.SetDelete(&delete_trexcLcLTTRExHVCluster);
      instance.SetDeleteArray(&deleteArray_trexcLcLTTRExHVCluster);
      instance.SetDestructor(&destruct_trexcLcLTTRExHVCluster);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::trex::TTRExHVCluster*)
   {
      return GenerateInitInstanceLocal((::trex::TTRExHVCluster*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::trex::TTRExHVCluster*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *trexcLcLTTRExHVCluster_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::trex::TTRExHVCluster*)0x0)->GetClass();
      trexcLcLTTRExHVCluster_TClassManip(theClass);
   return theClass;
   }

   static void trexcLcLTTRExHVCluster_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *trexcLcLTTRExPath_Dictionary();
   static void trexcLcLTTRExPath_TClassManip(TClass*);
   static void *new_trexcLcLTTRExPath(void *p = 0);
   static void *newArray_trexcLcLTTRExPath(Long_t size, void *p);
   static void delete_trexcLcLTTRExPath(void *p);
   static void deleteArray_trexcLcLTTRExPath(void *p);
   static void destruct_trexcLcLTTRExPath(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::trex::TTRExPath*)
   {
      ::trex::TTRExPath *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::trex::TTRExPath));
      static ::ROOT::TGenericClassInfo 
         instance("trex::TTRExPath", "TTRExPath.hxx", 41,
                  typeid(::trex::TTRExPath), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &trexcLcLTTRExPath_Dictionary, isa_proxy, 4,
                  sizeof(::trex::TTRExPath) );
      instance.SetNew(&new_trexcLcLTTRExPath);
      instance.SetNewArray(&newArray_trexcLcLTTRExPath);
      instance.SetDelete(&delete_trexcLcLTTRExPath);
      instance.SetDeleteArray(&deleteArray_trexcLcLTTRExPath);
      instance.SetDestructor(&destruct_trexcLcLTTRExPath);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::trex::TTRExPath*)
   {
      return GenerateInitInstanceLocal((::trex::TTRExPath*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::trex::TTRExPath*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *trexcLcLTTRExPath_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::trex::TTRExPath*)0x0)->GetClass();
      trexcLcLTTRExPath_TClassManip(theClass);
   return theClass;
   }

   static void trexcLcLTTRExPath_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_trexcLcLTTRExEvent(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::trex::TTRExEvent : new ::trex::TTRExEvent;
   }
   static void *newArray_trexcLcLTTRExEvent(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::trex::TTRExEvent[nElements] : new ::trex::TTRExEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_trexcLcLTTRExEvent(void *p) {
      delete ((::trex::TTRExEvent*)p);
   }
   static void deleteArray_trexcLcLTTRExEvent(void *p) {
      delete [] ((::trex::TTRExEvent*)p);
   }
   static void destruct_trexcLcLTTRExEvent(void *p) {
      typedef ::trex::TTRExEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::trex::TTRExEvent

namespace ROOT {
   // Wrappers around operator new
   static void *new_trexcLcLTTRExPattern(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::trex::TTRExPattern : new ::trex::TTRExPattern;
   }
   static void *newArray_trexcLcLTTRExPattern(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::trex::TTRExPattern[nElements] : new ::trex::TTRExPattern[nElements];
   }
   // Wrapper around operator delete
   static void delete_trexcLcLTTRExPattern(void *p) {
      delete ((::trex::TTRExPattern*)p);
   }
   static void deleteArray_trexcLcLTTRExPattern(void *p) {
      delete [] ((::trex::TTRExPattern*)p);
   }
   static void destruct_trexcLcLTTRExPattern(void *p) {
      typedef ::trex::TTRExPattern current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::trex::TTRExPattern

namespace ROOT {
   // Wrappers around operator new
   static void *new_trexcLcLTTPCHitPad(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::trex::TTPCHitPad : new ::trex::TTPCHitPad;
   }
   static void *newArray_trexcLcLTTPCHitPad(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::trex::TTPCHitPad[nElements] : new ::trex::TTPCHitPad[nElements];
   }
   // Wrapper around operator delete
   static void delete_trexcLcLTTPCHitPad(void *p) {
      delete ((::trex::TTPCHitPad*)p);
   }
   static void deleteArray_trexcLcLTTPCHitPad(void *p) {
      delete [] ((::trex::TTPCHitPad*)p);
   }
   static void destruct_trexcLcLTTPCHitPad(void *p) {
      typedef ::trex::TTPCHitPad current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::trex::TTPCHitPad

namespace ROOT {
   // Wrappers around operator new
   static void *new_trexcLcLTTRExHVCluster(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::trex::TTRExHVCluster : new ::trex::TTRExHVCluster;
   }
   static void *newArray_trexcLcLTTRExHVCluster(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::trex::TTRExHVCluster[nElements] : new ::trex::TTRExHVCluster[nElements];
   }
   // Wrapper around operator delete
   static void delete_trexcLcLTTRExHVCluster(void *p) {
      delete ((::trex::TTRExHVCluster*)p);
   }
   static void deleteArray_trexcLcLTTRExHVCluster(void *p) {
      delete [] ((::trex::TTRExHVCluster*)p);
   }
   static void destruct_trexcLcLTTRExHVCluster(void *p) {
      typedef ::trex::TTRExHVCluster current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::trex::TTRExHVCluster

namespace ROOT {
   // Wrappers around operator new
   static void *new_trexcLcLTTRExPath(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::trex::TTRExPath : new ::trex::TTRExPath;
   }
   static void *newArray_trexcLcLTTRExPath(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::trex::TTRExPath[nElements] : new ::trex::TTRExPath[nElements];
   }
   // Wrapper around operator delete
   static void delete_trexcLcLTTRExPath(void *p) {
      delete ((::trex::TTRExPath*)p);
   }
   static void deleteArray_trexcLcLTTRExPath(void *p) {
      delete [] ((::trex::TTRExPath*)p);
   }
   static void destruct_trexcLcLTTRExPath(void *p) {
      typedef ::trex::TTRExPath current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::trex::TTRExPath

namespace ROOT {
   static TClass *vectorlEvectorlEunsignedsPintgRsPgR_Dictionary();
   static void vectorlEvectorlEunsignedsPintgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEunsignedsPintgRsPgR(void *p = 0);
   static void *newArray_vectorlEvectorlEunsignedsPintgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEunsignedsPintgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEunsignedsPintgRsPgR(void *p);
   static void destruct_vectorlEvectorlEunsignedsPintgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<unsigned int> >*)
   {
      vector<vector<unsigned int> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<unsigned int> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<unsigned int> >", -2, "vector", 214,
                  typeid(vector<vector<unsigned int> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEunsignedsPintgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<vector<unsigned int> >) );
      instance.SetNew(&new_vectorlEvectorlEunsignedsPintgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEunsignedsPintgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEunsignedsPintgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEunsignedsPintgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEunsignedsPintgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<unsigned int> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<vector<unsigned int> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEunsignedsPintgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<unsigned int> >*)0x0)->GetClass();
      vectorlEvectorlEunsignedsPintgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEunsignedsPintgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEunsignedsPintgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<unsigned int> > : new vector<vector<unsigned int> >;
   }
   static void *newArray_vectorlEvectorlEunsignedsPintgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<unsigned int> >[nElements] : new vector<vector<unsigned int> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEunsignedsPintgRsPgR(void *p) {
      delete ((vector<vector<unsigned int> >*)p);
   }
   static void deleteArray_vectorlEvectorlEunsignedsPintgRsPgR(void *p) {
      delete [] ((vector<vector<unsigned int> >*)p);
   }
   static void destruct_vectorlEvectorlEunsignedsPintgRsPgR(void *p) {
      typedef vector<vector<unsigned int> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<unsigned int> >

namespace ROOT {
   static TClass *vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR_Dictionary();
   static void vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR(void *p = 0);
   static void *newArray_vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR(void *p);
   static void destruct_vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<trex::TTRExHVCluster> >*)
   {
      vector<vector<trex::TTRExHVCluster> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<trex::TTRExHVCluster> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<trex::TTRExHVCluster> >", -2, "vector", 214,
                  typeid(vector<vector<trex::TTRExHVCluster> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<trex::TTRExHVCluster> >) );
      instance.SetNew(&new_vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<trex::TTRExHVCluster> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<vector<trex::TTRExHVCluster> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<trex::TTRExHVCluster> >*)0x0)->GetClass();
      vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<trex::TTRExHVCluster> > : new vector<vector<trex::TTRExHVCluster> >;
   }
   static void *newArray_vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<trex::TTRExHVCluster> >[nElements] : new vector<vector<trex::TTRExHVCluster> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR(void *p) {
      delete ((vector<vector<trex::TTRExHVCluster> >*)p);
   }
   static void deleteArray_vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR(void *p) {
      delete [] ((vector<vector<trex::TTRExHVCluster> >*)p);
   }
   static void destruct_vectorlEvectorlEtrexcLcLTTRExHVClustergRsPgR(void *p) {
      typedef vector<vector<trex::TTRExHVCluster> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<trex::TTRExHVCluster> >

namespace ROOT {
   static TClass *vectorlEunsignedsPintgR_Dictionary();
   static void vectorlEunsignedsPintgR_TClassManip(TClass*);
   static void *new_vectorlEunsignedsPintgR(void *p = 0);
   static void *newArray_vectorlEunsignedsPintgR(Long_t size, void *p);
   static void delete_vectorlEunsignedsPintgR(void *p);
   static void deleteArray_vectorlEunsignedsPintgR(void *p);
   static void destruct_vectorlEunsignedsPintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<unsigned int>*)
   {
      vector<unsigned int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<unsigned int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<unsigned int>", -2, "vector", 214,
                  typeid(vector<unsigned int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEunsignedsPintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<unsigned int>) );
      instance.SetNew(&new_vectorlEunsignedsPintgR);
      instance.SetNewArray(&newArray_vectorlEunsignedsPintgR);
      instance.SetDelete(&delete_vectorlEunsignedsPintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEunsignedsPintgR);
      instance.SetDestructor(&destruct_vectorlEunsignedsPintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<unsigned int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<unsigned int>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEunsignedsPintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<unsigned int>*)0x0)->GetClass();
      vectorlEunsignedsPintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEunsignedsPintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEunsignedsPintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned int> : new vector<unsigned int>;
   }
   static void *newArray_vectorlEunsignedsPintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned int>[nElements] : new vector<unsigned int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEunsignedsPintgR(void *p) {
      delete ((vector<unsigned int>*)p);
   }
   static void deleteArray_vectorlEunsignedsPintgR(void *p) {
      delete [] ((vector<unsigned int>*)p);
   }
   static void destruct_vectorlEunsignedsPintgR(void *p) {
      typedef vector<unsigned int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<unsigned int>

namespace ROOT {
   static TClass *vectorlEtrexcLcLTTRExPatterngR_Dictionary();
   static void vectorlEtrexcLcLTTRExPatterngR_TClassManip(TClass*);
   static void *new_vectorlEtrexcLcLTTRExPatterngR(void *p = 0);
   static void *newArray_vectorlEtrexcLcLTTRExPatterngR(Long_t size, void *p);
   static void delete_vectorlEtrexcLcLTTRExPatterngR(void *p);
   static void deleteArray_vectorlEtrexcLcLTTRExPatterngR(void *p);
   static void destruct_vectorlEtrexcLcLTTRExPatterngR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<trex::TTRExPattern>*)
   {
      vector<trex::TTRExPattern> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<trex::TTRExPattern>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<trex::TTRExPattern>", -2, "vector", 214,
                  typeid(vector<trex::TTRExPattern>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEtrexcLcLTTRExPatterngR_Dictionary, isa_proxy, 4,
                  sizeof(vector<trex::TTRExPattern>) );
      instance.SetNew(&new_vectorlEtrexcLcLTTRExPatterngR);
      instance.SetNewArray(&newArray_vectorlEtrexcLcLTTRExPatterngR);
      instance.SetDelete(&delete_vectorlEtrexcLcLTTRExPatterngR);
      instance.SetDeleteArray(&deleteArray_vectorlEtrexcLcLTTRExPatterngR);
      instance.SetDestructor(&destruct_vectorlEtrexcLcLTTRExPatterngR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<trex::TTRExPattern> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<trex::TTRExPattern>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEtrexcLcLTTRExPatterngR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<trex::TTRExPattern>*)0x0)->GetClass();
      vectorlEtrexcLcLTTRExPatterngR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEtrexcLcLTTRExPatterngR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEtrexcLcLTTRExPatterngR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<trex::TTRExPattern> : new vector<trex::TTRExPattern>;
   }
   static void *newArray_vectorlEtrexcLcLTTRExPatterngR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<trex::TTRExPattern>[nElements] : new vector<trex::TTRExPattern>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEtrexcLcLTTRExPatterngR(void *p) {
      delete ((vector<trex::TTRExPattern>*)p);
   }
   static void deleteArray_vectorlEtrexcLcLTTRExPatterngR(void *p) {
      delete [] ((vector<trex::TTRExPattern>*)p);
   }
   static void destruct_vectorlEtrexcLcLTTRExPatterngR(void *p) {
      typedef vector<trex::TTRExPattern> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<trex::TTRExPattern>

namespace ROOT {
   static TClass *vectorlEtrexcLcLTTRExPathgR_Dictionary();
   static void vectorlEtrexcLcLTTRExPathgR_TClassManip(TClass*);
   static void *new_vectorlEtrexcLcLTTRExPathgR(void *p = 0);
   static void *newArray_vectorlEtrexcLcLTTRExPathgR(Long_t size, void *p);
   static void delete_vectorlEtrexcLcLTTRExPathgR(void *p);
   static void deleteArray_vectorlEtrexcLcLTTRExPathgR(void *p);
   static void destruct_vectorlEtrexcLcLTTRExPathgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<trex::TTRExPath>*)
   {
      vector<trex::TTRExPath> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<trex::TTRExPath>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<trex::TTRExPath>", -2, "vector", 214,
                  typeid(vector<trex::TTRExPath>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEtrexcLcLTTRExPathgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<trex::TTRExPath>) );
      instance.SetNew(&new_vectorlEtrexcLcLTTRExPathgR);
      instance.SetNewArray(&newArray_vectorlEtrexcLcLTTRExPathgR);
      instance.SetDelete(&delete_vectorlEtrexcLcLTTRExPathgR);
      instance.SetDeleteArray(&deleteArray_vectorlEtrexcLcLTTRExPathgR);
      instance.SetDestructor(&destruct_vectorlEtrexcLcLTTRExPathgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<trex::TTRExPath> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<trex::TTRExPath>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEtrexcLcLTTRExPathgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<trex::TTRExPath>*)0x0)->GetClass();
      vectorlEtrexcLcLTTRExPathgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEtrexcLcLTTRExPathgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEtrexcLcLTTRExPathgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<trex::TTRExPath> : new vector<trex::TTRExPath>;
   }
   static void *newArray_vectorlEtrexcLcLTTRExPathgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<trex::TTRExPath>[nElements] : new vector<trex::TTRExPath>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEtrexcLcLTTRExPathgR(void *p) {
      delete ((vector<trex::TTRExPath>*)p);
   }
   static void deleteArray_vectorlEtrexcLcLTTRExPathgR(void *p) {
      delete [] ((vector<trex::TTRExPath>*)p);
   }
   static void destruct_vectorlEtrexcLcLTTRExPathgR(void *p) {
      typedef vector<trex::TTRExPath> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<trex::TTRExPath>

namespace ROOT {
   static TClass *vectorlEtrexcLcLTTRExJunctiongR_Dictionary();
   static void vectorlEtrexcLcLTTRExJunctiongR_TClassManip(TClass*);
   static void *new_vectorlEtrexcLcLTTRExJunctiongR(void *p = 0);
   static void *newArray_vectorlEtrexcLcLTTRExJunctiongR(Long_t size, void *p);
   static void delete_vectorlEtrexcLcLTTRExJunctiongR(void *p);
   static void deleteArray_vectorlEtrexcLcLTTRExJunctiongR(void *p);
   static void destruct_vectorlEtrexcLcLTTRExJunctiongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<trex::TTRExJunction>*)
   {
      vector<trex::TTRExJunction> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<trex::TTRExJunction>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<trex::TTRExJunction>", -2, "vector", 214,
                  typeid(vector<trex::TTRExJunction>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEtrexcLcLTTRExJunctiongR_Dictionary, isa_proxy, 0,
                  sizeof(vector<trex::TTRExJunction>) );
      instance.SetNew(&new_vectorlEtrexcLcLTTRExJunctiongR);
      instance.SetNewArray(&newArray_vectorlEtrexcLcLTTRExJunctiongR);
      instance.SetDelete(&delete_vectorlEtrexcLcLTTRExJunctiongR);
      instance.SetDeleteArray(&deleteArray_vectorlEtrexcLcLTTRExJunctiongR);
      instance.SetDestructor(&destruct_vectorlEtrexcLcLTTRExJunctiongR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<trex::TTRExJunction> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<trex::TTRExJunction>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEtrexcLcLTTRExJunctiongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<trex::TTRExJunction>*)0x0)->GetClass();
      vectorlEtrexcLcLTTRExJunctiongR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEtrexcLcLTTRExJunctiongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEtrexcLcLTTRExJunctiongR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<trex::TTRExJunction> : new vector<trex::TTRExJunction>;
   }
   static void *newArray_vectorlEtrexcLcLTTRExJunctiongR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<trex::TTRExJunction>[nElements] : new vector<trex::TTRExJunction>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEtrexcLcLTTRExJunctiongR(void *p) {
      delete ((vector<trex::TTRExJunction>*)p);
   }
   static void deleteArray_vectorlEtrexcLcLTTRExJunctiongR(void *p) {
      delete [] ((vector<trex::TTRExJunction>*)p);
   }
   static void destruct_vectorlEtrexcLcLTTRExJunctiongR(void *p) {
      typedef vector<trex::TTRExJunction> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<trex::TTRExJunction>

namespace ROOT {
   static TClass *vectorlEtrexcLcLTTRExHVClustergR_Dictionary();
   static void vectorlEtrexcLcLTTRExHVClustergR_TClassManip(TClass*);
   static void *new_vectorlEtrexcLcLTTRExHVClustergR(void *p = 0);
   static void *newArray_vectorlEtrexcLcLTTRExHVClustergR(Long_t size, void *p);
   static void delete_vectorlEtrexcLcLTTRExHVClustergR(void *p);
   static void deleteArray_vectorlEtrexcLcLTTRExHVClustergR(void *p);
   static void destruct_vectorlEtrexcLcLTTRExHVClustergR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<trex::TTRExHVCluster>*)
   {
      vector<trex::TTRExHVCluster> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<trex::TTRExHVCluster>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<trex::TTRExHVCluster>", -2, "vector", 214,
                  typeid(vector<trex::TTRExHVCluster>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEtrexcLcLTTRExHVClustergR_Dictionary, isa_proxy, 4,
                  sizeof(vector<trex::TTRExHVCluster>) );
      instance.SetNew(&new_vectorlEtrexcLcLTTRExHVClustergR);
      instance.SetNewArray(&newArray_vectorlEtrexcLcLTTRExHVClustergR);
      instance.SetDelete(&delete_vectorlEtrexcLcLTTRExHVClustergR);
      instance.SetDeleteArray(&deleteArray_vectorlEtrexcLcLTTRExHVClustergR);
      instance.SetDestructor(&destruct_vectorlEtrexcLcLTTRExHVClustergR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<trex::TTRExHVCluster> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<trex::TTRExHVCluster>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEtrexcLcLTTRExHVClustergR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<trex::TTRExHVCluster>*)0x0)->GetClass();
      vectorlEtrexcLcLTTRExHVClustergR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEtrexcLcLTTRExHVClustergR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEtrexcLcLTTRExHVClustergR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<trex::TTRExHVCluster> : new vector<trex::TTRExHVCluster>;
   }
   static void *newArray_vectorlEtrexcLcLTTRExHVClustergR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<trex::TTRExHVCluster>[nElements] : new vector<trex::TTRExHVCluster>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEtrexcLcLTTRExHVClustergR(void *p) {
      delete ((vector<trex::TTRExHVCluster>*)p);
   }
   static void deleteArray_vectorlEtrexcLcLTTRExHVClustergR(void *p) {
      delete [] ((vector<trex::TTRExHVCluster>*)p);
   }
   static void destruct_vectorlEtrexcLcLTTRExHVClustergR(void *p) {
      typedef vector<trex::TTRExHVCluster> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<trex::TTRExHVCluster>

namespace ROOT {
   static TClass *vectorlEtrexcLcLTTPCHitPadgR_Dictionary();
   static void vectorlEtrexcLcLTTPCHitPadgR_TClassManip(TClass*);
   static void *new_vectorlEtrexcLcLTTPCHitPadgR(void *p = 0);
   static void *newArray_vectorlEtrexcLcLTTPCHitPadgR(Long_t size, void *p);
   static void delete_vectorlEtrexcLcLTTPCHitPadgR(void *p);
   static void deleteArray_vectorlEtrexcLcLTTPCHitPadgR(void *p);
   static void destruct_vectorlEtrexcLcLTTPCHitPadgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<trex::TTPCHitPad>*)
   {
      vector<trex::TTPCHitPad> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<trex::TTPCHitPad>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<trex::TTPCHitPad>", -2, "vector", 214,
                  typeid(vector<trex::TTPCHitPad>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEtrexcLcLTTPCHitPadgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<trex::TTPCHitPad>) );
      instance.SetNew(&new_vectorlEtrexcLcLTTPCHitPadgR);
      instance.SetNewArray(&newArray_vectorlEtrexcLcLTTPCHitPadgR);
      instance.SetDelete(&delete_vectorlEtrexcLcLTTPCHitPadgR);
      instance.SetDeleteArray(&deleteArray_vectorlEtrexcLcLTTPCHitPadgR);
      instance.SetDestructor(&destruct_vectorlEtrexcLcLTTPCHitPadgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<trex::TTPCHitPad> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<trex::TTPCHitPad>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEtrexcLcLTTPCHitPadgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<trex::TTPCHitPad>*)0x0)->GetClass();
      vectorlEtrexcLcLTTPCHitPadgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEtrexcLcLTTPCHitPadgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEtrexcLcLTTPCHitPadgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<trex::TTPCHitPad> : new vector<trex::TTPCHitPad>;
   }
   static void *newArray_vectorlEtrexcLcLTTPCHitPadgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<trex::TTPCHitPad>[nElements] : new vector<trex::TTPCHitPad>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEtrexcLcLTTPCHitPadgR(void *p) {
      delete ((vector<trex::TTPCHitPad>*)p);
   }
   static void deleteArray_vectorlEtrexcLcLTTPCHitPadgR(void *p) {
      delete [] ((vector<trex::TTPCHitPad>*)p);
   }
   static void destruct_vectorlEtrexcLcLTTPCHitPadgR(void *p) {
      typedef vector<trex::TTPCHitPad> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<trex::TTPCHitPad>

namespace ROOT {
   static TClass *vectorlEtrexcLcLTTPCHitPadmUgR_Dictionary();
   static void vectorlEtrexcLcLTTPCHitPadmUgR_TClassManip(TClass*);
   static void *new_vectorlEtrexcLcLTTPCHitPadmUgR(void *p = 0);
   static void *newArray_vectorlEtrexcLcLTTPCHitPadmUgR(Long_t size, void *p);
   static void delete_vectorlEtrexcLcLTTPCHitPadmUgR(void *p);
   static void deleteArray_vectorlEtrexcLcLTTPCHitPadmUgR(void *p);
   static void destruct_vectorlEtrexcLcLTTPCHitPadmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<trex::TTPCHitPad*>*)
   {
      vector<trex::TTPCHitPad*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<trex::TTPCHitPad*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<trex::TTPCHitPad*>", -2, "vector", 214,
                  typeid(vector<trex::TTPCHitPad*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEtrexcLcLTTPCHitPadmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<trex::TTPCHitPad*>) );
      instance.SetNew(&new_vectorlEtrexcLcLTTPCHitPadmUgR);
      instance.SetNewArray(&newArray_vectorlEtrexcLcLTTPCHitPadmUgR);
      instance.SetDelete(&delete_vectorlEtrexcLcLTTPCHitPadmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEtrexcLcLTTPCHitPadmUgR);
      instance.SetDestructor(&destruct_vectorlEtrexcLcLTTPCHitPadmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<trex::TTPCHitPad*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<trex::TTPCHitPad*>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEtrexcLcLTTPCHitPadmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<trex::TTPCHitPad*>*)0x0)->GetClass();
      vectorlEtrexcLcLTTPCHitPadmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEtrexcLcLTTPCHitPadmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEtrexcLcLTTPCHitPadmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<trex::TTPCHitPad*> : new vector<trex::TTPCHitPad*>;
   }
   static void *newArray_vectorlEtrexcLcLTTPCHitPadmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<trex::TTPCHitPad*>[nElements] : new vector<trex::TTPCHitPad*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEtrexcLcLTTPCHitPadmUgR(void *p) {
      delete ((vector<trex::TTPCHitPad*>*)p);
   }
   static void deleteArray_vectorlEtrexcLcLTTPCHitPadmUgR(void *p) {
      delete [] ((vector<trex::TTPCHitPad*>*)p);
   }
   static void destruct_vectorlEtrexcLcLTTPCHitPadmUgR(void *p) {
      typedef vector<trex::TTPCHitPad*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<trex::TTPCHitPad*>

namespace {
  void TriggerDictionaryInitialization_TTRExPatternDict_Impl() {
    static const char* headers[] = {
"TTRExPattern.hxx",
"OutLink.hh",
0
    };
    static const char* includePaths[] = {
"/storage/epp2/epp_software/2015.1/Cellar/root6/6.06.06/include/root",
"/home/epp/phsmaj/hptrex/src/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TTRExPatternDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace trex{class __attribute__((annotate("$clingAutoload$TTRExPattern.hxx")))  TTRExEvent;}
namespace trex{class __attribute__((annotate("$clingAutoload$TTRExPattern.hxx")))  TTRExPattern;}
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$string")))  allocator;
}
namespace trex{class __attribute__((annotate("$clingAutoload$TTRExPattern.hxx")))  TTPCHitPad;}
namespace trex{class __attribute__((annotate("$clingAutoload$TTRExPattern.hxx")))  TTRExHVCluster;}
namespace trex{class __attribute__((annotate("$clingAutoload$TTRExPattern.hxx")))  TTRExPath;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TTRExPatternDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TTRExPattern.hxx"
#include "OutLink.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"trex::TTPCHitPad", payloadCode, "@",
"trex::TTRExEvent", payloadCode, "@",
"trex::TTRExHVCluster", payloadCode, "@",
"trex::TTRExPath", payloadCode, "@",
"trex::TTRExPattern", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TTRExPatternDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TTRExPatternDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TTRExPatternDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TTRExPatternDict() {
  TriggerDictionaryInitialization_TTRExPatternDict_Impl();
}
