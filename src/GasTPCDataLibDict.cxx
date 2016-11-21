// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME GasTPCDataLibDict

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
#include "GasTPCDataLib.hxx"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_namedRecord(void *p = 0);
   static void *newArray_namedRecord(Long_t size, void *p);
   static void delete_namedRecord(void *p);
   static void deleteArray_namedRecord(void *p);
   static void destruct_namedRecord(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::namedRecord*)
   {
      ::namedRecord *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::namedRecord >(0);
      static ::ROOT::TGenericClassInfo 
         instance("namedRecord", ::namedRecord::Class_Version(), "GasTPCDataLib.hxx", 34,
                  typeid(::namedRecord), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::namedRecord::Dictionary, isa_proxy, 4,
                  sizeof(::namedRecord) );
      instance.SetNew(&new_namedRecord);
      instance.SetNewArray(&newArray_namedRecord);
      instance.SetDelete(&delete_namedRecord);
      instance.SetDeleteArray(&deleteArray_namedRecord);
      instance.SetDestructor(&destruct_namedRecord);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::namedRecord*)
   {
      return GenerateInitInstanceLocal((::namedRecord*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::namedRecord*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ParticleDescrShortRecord(void *p = 0);
   static void *newArray_ParticleDescrShortRecord(Long_t size, void *p);
   static void delete_ParticleDescrShortRecord(void *p);
   static void deleteArray_ParticleDescrShortRecord(void *p);
   static void destruct_ParticleDescrShortRecord(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ParticleDescrShortRecord*)
   {
      ::ParticleDescrShortRecord *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ParticleDescrShortRecord >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ParticleDescrShortRecord", ::ParticleDescrShortRecord::Class_Version(), "GasTPCDataLib.hxx", 48,
                  typeid(::ParticleDescrShortRecord), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ParticleDescrShortRecord::Dictionary, isa_proxy, 4,
                  sizeof(::ParticleDescrShortRecord) );
      instance.SetNew(&new_ParticleDescrShortRecord);
      instance.SetNewArray(&newArray_ParticleDescrShortRecord);
      instance.SetDelete(&delete_ParticleDescrShortRecord);
      instance.SetDeleteArray(&deleteArray_ParticleDescrShortRecord);
      instance.SetDestructor(&destruct_ParticleDescrShortRecord);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ParticleDescrShortRecord*)
   {
      return GenerateInitInstanceLocal((::ParticleDescrShortRecord*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ParticleDescrShortRecord*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ParticleDescrRecord(void *p = 0);
   static void *newArray_ParticleDescrRecord(Long_t size, void *p);
   static void delete_ParticleDescrRecord(void *p);
   static void deleteArray_ParticleDescrRecord(void *p);
   static void destruct_ParticleDescrRecord(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ParticleDescrRecord*)
   {
      ::ParticleDescrRecord *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ParticleDescrRecord >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ParticleDescrRecord", ::ParticleDescrRecord::Class_Version(), "GasTPCDataLib.hxx", 80,
                  typeid(::ParticleDescrRecord), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ParticleDescrRecord::Dictionary, isa_proxy, 4,
                  sizeof(::ParticleDescrRecord) );
      instance.SetNew(&new_ParticleDescrRecord);
      instance.SetNewArray(&newArray_ParticleDescrRecord);
      instance.SetDelete(&delete_ParticleDescrRecord);
      instance.SetDeleteArray(&deleteArray_ParticleDescrRecord);
      instance.SetDestructor(&destruct_ParticleDescrRecord);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ParticleDescrRecord*)
   {
      return GenerateInitInstanceLocal((::ParticleDescrRecord*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ParticleDescrRecord*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_BaseEventRecord(void *p = 0);
   static void *newArray_BaseEventRecord(Long_t size, void *p);
   static void delete_BaseEventRecord(void *p);
   static void deleteArray_BaseEventRecord(void *p);
   static void destruct_BaseEventRecord(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::BaseEventRecord*)
   {
      ::BaseEventRecord *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::BaseEventRecord >(0);
      static ::ROOT::TGenericClassInfo 
         instance("BaseEventRecord", ::BaseEventRecord::Class_Version(), "GasTPCDataLib.hxx", 103,
                  typeid(::BaseEventRecord), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::BaseEventRecord::Dictionary, isa_proxy, 4,
                  sizeof(::BaseEventRecord) );
      instance.SetNew(&new_BaseEventRecord);
      instance.SetNewArray(&newArray_BaseEventRecord);
      instance.SetDelete(&delete_BaseEventRecord);
      instance.SetDeleteArray(&deleteArray_BaseEventRecord);
      instance.SetDestructor(&destruct_BaseEventRecord);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::BaseEventRecord*)
   {
      return GenerateInitInstanceLocal((::BaseEventRecord*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::BaseEventRecord*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_NeutrinoHit(void *p = 0);
   static void *newArray_NeutrinoHit(Long_t size, void *p);
   static void delete_NeutrinoHit(void *p);
   static void deleteArray_NeutrinoHit(void *p);
   static void destruct_NeutrinoHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NeutrinoHit*)
   {
      ::NeutrinoHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::NeutrinoHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("NeutrinoHit", ::NeutrinoHit::Class_Version(), "GasTPCDataLib.hxx", 134,
                  typeid(::NeutrinoHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::NeutrinoHit::Dictionary, isa_proxy, 4,
                  sizeof(::NeutrinoHit) );
      instance.SetNew(&new_NeutrinoHit);
      instance.SetNewArray(&newArray_NeutrinoHit);
      instance.SetDelete(&delete_NeutrinoHit);
      instance.SetDeleteArray(&deleteArray_NeutrinoHit);
      instance.SetDestructor(&destruct_NeutrinoHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NeutrinoHit*)
   {
      return GenerateInitInstanceLocal((::NeutrinoHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::NeutrinoHit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_MuonDecayEvent(void *p = 0);
   static void *newArray_MuonDecayEvent(Long_t size, void *p);
   static void delete_MuonDecayEvent(void *p);
   static void deleteArray_MuonDecayEvent(void *p);
   static void destruct_MuonDecayEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MuonDecayEvent*)
   {
      ::MuonDecayEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MuonDecayEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MuonDecayEvent", ::MuonDecayEvent::Class_Version(), "GasTPCDataLib.hxx", 192,
                  typeid(::MuonDecayEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MuonDecayEvent::Dictionary, isa_proxy, 4,
                  sizeof(::MuonDecayEvent) );
      instance.SetNew(&new_MuonDecayEvent);
      instance.SetNewArray(&newArray_MuonDecayEvent);
      instance.SetDelete(&delete_MuonDecayEvent);
      instance.SetDeleteArray(&deleteArray_MuonDecayEvent);
      instance.SetDestructor(&destruct_MuonDecayEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MuonDecayEvent*)
   {
      return GenerateInitInstanceLocal((::MuonDecayEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::MuonDecayEvent*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PionDecayEvent(void *p = 0);
   static void *newArray_PionDecayEvent(Long_t size, void *p);
   static void delete_PionDecayEvent(void *p);
   static void deleteArray_PionDecayEvent(void *p);
   static void destruct_PionDecayEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PionDecayEvent*)
   {
      ::PionDecayEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PionDecayEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PionDecayEvent", ::PionDecayEvent::Class_Version(), "GasTPCDataLib.hxx", 157,
                  typeid(::PionDecayEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PionDecayEvent::Dictionary, isa_proxy, 4,
                  sizeof(::PionDecayEvent) );
      instance.SetNew(&new_PionDecayEvent);
      instance.SetNewArray(&newArray_PionDecayEvent);
      instance.SetDelete(&delete_PionDecayEvent);
      instance.SetDeleteArray(&deleteArray_PionDecayEvent);
      instance.SetDestructor(&destruct_PionDecayEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PionDecayEvent*)
   {
      return GenerateInitInstanceLocal((::PionDecayEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::PionDecayEvent*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_NeutrinoEvent(void *p = 0);
   static void *newArray_NeutrinoEvent(Long_t size, void *p);
   static void delete_NeutrinoEvent(void *p);
   static void deleteArray_NeutrinoEvent(void *p);
   static void destruct_NeutrinoEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NeutrinoEvent*)
   {
      ::NeutrinoEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::NeutrinoEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("NeutrinoEvent", ::NeutrinoEvent::Class_Version(), "GasTPCDataLib.hxx", 228,
                  typeid(::NeutrinoEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::NeutrinoEvent::Dictionary, isa_proxy, 4,
                  sizeof(::NeutrinoEvent) );
      instance.SetNew(&new_NeutrinoEvent);
      instance.SetNewArray(&newArray_NeutrinoEvent);
      instance.SetDelete(&delete_NeutrinoEvent);
      instance.SetDeleteArray(&deleteArray_NeutrinoEvent);
      instance.SetDestructor(&destruct_NeutrinoEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NeutrinoEvent*)
   {
      return GenerateInitInstanceLocal((::NeutrinoEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::NeutrinoEvent*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SDHit(void *p = 0);
   static void *newArray_SDHit(Long_t size, void *p);
   static void delete_SDHit(void *p);
   static void deleteArray_SDHit(void *p);
   static void destruct_SDHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SDHit*)
   {
      ::SDHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SDHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SDHit", ::SDHit::Class_Version(), "GasTPCDataLib.hxx", 531,
                  typeid(::SDHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SDHit::Dictionary, isa_proxy, 4,
                  sizeof(::SDHit) );
      instance.SetNew(&new_SDHit);
      instance.SetNewArray(&newArray_SDHit);
      instance.SetDelete(&delete_SDHit);
      instance.SetDeleteArray(&deleteArray_SDHit);
      instance.SetDestructor(&destruct_SDHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SDHit*)
   {
      return GenerateInitInstanceLocal((::SDHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::SDHit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_tpcFidHit(void *p = 0);
   static void *newArray_tpcFidHit(Long_t size, void *p);
   static void delete_tpcFidHit(void *p);
   static void deleteArray_tpcFidHit(void *p);
   static void destruct_tpcFidHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::tpcFidHit*)
   {
      ::tpcFidHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::tpcFidHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("tpcFidHit", ::tpcFidHit::Class_Version(), "GasTPCDataLib.hxx", 618,
                  typeid(::tpcFidHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::tpcFidHit::Dictionary, isa_proxy, 4,
                  sizeof(::tpcFidHit) );
      instance.SetNew(&new_tpcFidHit);
      instance.SetNewArray(&newArray_tpcFidHit);
      instance.SetDelete(&delete_tpcFidHit);
      instance.SetDeleteArray(&deleteArray_tpcFidHit);
      instance.SetDestructor(&destruct_tpcFidHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::tpcFidHit*)
   {
      return GenerateInitInstanceLocal((::tpcFidHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::tpcFidHit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_scintHit(void *p = 0);
   static void *newArray_scintHit(Long_t size, void *p);
   static void delete_scintHit(void *p);
   static void deleteArray_scintHit(void *p);
   static void destruct_scintHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::scintHit*)
   {
      ::scintHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::scintHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("scintHit", ::scintHit::Class_Version(), "GasTPCDataLib.hxx", 631,
                  typeid(::scintHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::scintHit::Dictionary, isa_proxy, 4,
                  sizeof(::scintHit) );
      instance.SetNew(&new_scintHit);
      instance.SetNewArray(&newArray_scintHit);
      instance.SetDelete(&delete_scintHit);
      instance.SetDeleteArray(&deleteArray_scintHit);
      instance.SetDestructor(&destruct_scintHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::scintHit*)
   {
      return GenerateInitInstanceLocal((::scintHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::scintHit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_MINDHit(void *p = 0);
   static void *newArray_MINDHit(Long_t size, void *p);
   static void delete_MINDHit(void *p);
   static void deleteArray_MINDHit(void *p);
   static void destruct_MINDHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MINDHit*)
   {
      ::MINDHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MINDHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MINDHit", ::MINDHit::Class_Version(), "GasTPCDataLib.hxx", 654,
                  typeid(::MINDHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MINDHit::Dictionary, isa_proxy, 4,
                  sizeof(::MINDHit) );
      instance.SetNew(&new_MINDHit);
      instance.SetNewArray(&newArray_MINDHit);
      instance.SetDelete(&delete_MINDHit);
      instance.SetDeleteArray(&deleteArray_MINDHit);
      instance.SetDestructor(&destruct_MINDHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MINDHit*)
   {
      return GenerateInitInstanceLocal((::MINDHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::MINDHit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SimulData(void *p = 0);
   static void *newArray_SimulData(Long_t size, void *p);
   static void delete_SimulData(void *p);
   static void deleteArray_SimulData(void *p);
   static void destruct_SimulData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SimulData*)
   {
      ::SimulData *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SimulData >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SimulData", ::SimulData::Class_Version(), "GasTPCDataLib.hxx", 674,
                  typeid(::SimulData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SimulData::Dictionary, isa_proxy, 4,
                  sizeof(::SimulData) );
      instance.SetNew(&new_SimulData);
      instance.SetNewArray(&newArray_SimulData);
      instance.SetDelete(&delete_SimulData);
      instance.SetDeleteArray(&deleteArray_SimulData);
      instance.SetDestructor(&destruct_SimulData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SimulData*)
   {
      return GenerateInitInstanceLocal((::SimulData*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::SimulData*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_GeantBasicParticle(void *p = 0);
   static void *newArray_GeantBasicParticle(Long_t size, void *p);
   static void delete_GeantBasicParticle(void *p);
   static void deleteArray_GeantBasicParticle(void *p);
   static void destruct_GeantBasicParticle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::GeantBasicParticle*)
   {
      ::GeantBasicParticle *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::GeantBasicParticle >(0);
      static ::ROOT::TGenericClassInfo 
         instance("GeantBasicParticle", ::GeantBasicParticle::Class_Version(), "GasTPCDataLib.hxx", 313,
                  typeid(::GeantBasicParticle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::GeantBasicParticle::Dictionary, isa_proxy, 4,
                  sizeof(::GeantBasicParticle) );
      instance.SetNew(&new_GeantBasicParticle);
      instance.SetNewArray(&newArray_GeantBasicParticle);
      instance.SetDelete(&delete_GeantBasicParticle);
      instance.SetDeleteArray(&deleteArray_GeantBasicParticle);
      instance.SetDestructor(&destruct_GeantBasicParticle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::GeantBasicParticle*)
   {
      return GenerateInitInstanceLocal((::GeantBasicParticle*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::GeantBasicParticle*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_GeantParticle(void *p = 0);
   static void *newArray_GeantParticle(Long_t size, void *p);
   static void delete_GeantParticle(void *p);
   static void deleteArray_GeantParticle(void *p);
   static void destruct_GeantParticle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::GeantParticle*)
   {
      ::GeantParticle *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::GeantParticle >(0);
      static ::ROOT::TGenericClassInfo 
         instance("GeantParticle", ::GeantParticle::Class_Version(), "GasTPCDataLib.hxx", 358,
                  typeid(::GeantParticle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::GeantParticle::Dictionary, isa_proxy, 4,
                  sizeof(::GeantParticle) );
      instance.SetNew(&new_GeantParticle);
      instance.SetNewArray(&newArray_GeantParticle);
      instance.SetDelete(&delete_GeantParticle);
      instance.SetDeleteArray(&deleteArray_GeantParticle);
      instance.SetDestructor(&destruct_GeantParticle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::GeantParticle*)
   {
      return GenerateInitInstanceLocal((::GeantParticle*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::GeantParticle*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_GeantTrackingTruth(void *p = 0);
   static void *newArray_GeantTrackingTruth(Long_t size, void *p);
   static void delete_GeantTrackingTruth(void *p);
   static void deleteArray_GeantTrackingTruth(void *p);
   static void destruct_GeantTrackingTruth(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::GeantTrackingTruth*)
   {
      ::GeantTrackingTruth *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::GeantTrackingTruth >(0);
      static ::ROOT::TGenericClassInfo 
         instance("GeantTrackingTruth", ::GeantTrackingTruth::Class_Version(), "GasTPCDataLib.hxx", 483,
                  typeid(::GeantTrackingTruth), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::GeantTrackingTruth::Dictionary, isa_proxy, 4,
                  sizeof(::GeantTrackingTruth) );
      instance.SetNew(&new_GeantTrackingTruth);
      instance.SetNewArray(&newArray_GeantTrackingTruth);
      instance.SetDelete(&delete_GeantTrackingTruth);
      instance.SetDeleteArray(&deleteArray_GeantTrackingTruth);
      instance.SetDestructor(&destruct_GeantTrackingTruth);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::GeantTrackingTruth*)
   {
      return GenerateInitInstanceLocal((::GeantTrackingTruth*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::GeantTrackingTruth*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_MemTestRecord(void *p = 0);
   static void *newArray_MemTestRecord(Long_t size, void *p);
   static void delete_MemTestRecord(void *p);
   static void deleteArray_MemTestRecord(void *p);
   static void destruct_MemTestRecord(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MemTestRecord*)
   {
      ::MemTestRecord *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MemTestRecord >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MemTestRecord", ::MemTestRecord::Class_Version(), "GasTPCDataLib.hxx", 709,
                  typeid(::MemTestRecord), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MemTestRecord::Dictionary, isa_proxy, 4,
                  sizeof(::MemTestRecord) );
      instance.SetNew(&new_MemTestRecord);
      instance.SetNewArray(&newArray_MemTestRecord);
      instance.SetDelete(&delete_MemTestRecord);
      instance.SetDeleteArray(&deleteArray_MemTestRecord);
      instance.SetDestructor(&destruct_MemTestRecord);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MemTestRecord*)
   {
      return GenerateInitInstanceLocal((::MemTestRecord*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::MemTestRecord*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr namedRecord::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *namedRecord::Class_Name()
{
   return "namedRecord";
}

//______________________________________________________________________________
const char *namedRecord::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::namedRecord*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int namedRecord::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::namedRecord*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *namedRecord::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::namedRecord*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *namedRecord::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::namedRecord*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ParticleDescrShortRecord::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ParticleDescrShortRecord::Class_Name()
{
   return "ParticleDescrShortRecord";
}

//______________________________________________________________________________
const char *ParticleDescrShortRecord::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ParticleDescrShortRecord*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ParticleDescrShortRecord::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ParticleDescrShortRecord*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ParticleDescrShortRecord::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ParticleDescrShortRecord*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ParticleDescrShortRecord::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ParticleDescrShortRecord*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ParticleDescrRecord::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ParticleDescrRecord::Class_Name()
{
   return "ParticleDescrRecord";
}

//______________________________________________________________________________
const char *ParticleDescrRecord::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ParticleDescrRecord*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ParticleDescrRecord::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ParticleDescrRecord*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ParticleDescrRecord::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ParticleDescrRecord*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ParticleDescrRecord::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ParticleDescrRecord*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr BaseEventRecord::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *BaseEventRecord::Class_Name()
{
   return "BaseEventRecord";
}

//______________________________________________________________________________
const char *BaseEventRecord::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::BaseEventRecord*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int BaseEventRecord::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::BaseEventRecord*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *BaseEventRecord::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::BaseEventRecord*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *BaseEventRecord::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::BaseEventRecord*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr NeutrinoHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *NeutrinoHit::Class_Name()
{
   return "NeutrinoHit";
}

//______________________________________________________________________________
const char *NeutrinoHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NeutrinoHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int NeutrinoHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NeutrinoHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *NeutrinoHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NeutrinoHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *NeutrinoHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NeutrinoHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr MuonDecayEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MuonDecayEvent::Class_Name()
{
   return "MuonDecayEvent";
}

//______________________________________________________________________________
const char *MuonDecayEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MuonDecayEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MuonDecayEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MuonDecayEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MuonDecayEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MuonDecayEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MuonDecayEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MuonDecayEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr PionDecayEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PionDecayEvent::Class_Name()
{
   return "PionDecayEvent";
}

//______________________________________________________________________________
const char *PionDecayEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PionDecayEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PionDecayEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PionDecayEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PionDecayEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PionDecayEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PionDecayEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PionDecayEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr NeutrinoEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *NeutrinoEvent::Class_Name()
{
   return "NeutrinoEvent";
}

//______________________________________________________________________________
const char *NeutrinoEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NeutrinoEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int NeutrinoEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NeutrinoEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *NeutrinoEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NeutrinoEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *NeutrinoEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NeutrinoEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SDHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SDHit::Class_Name()
{
   return "SDHit";
}

//______________________________________________________________________________
const char *SDHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SDHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SDHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SDHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SDHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SDHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SDHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SDHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr tpcFidHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *tpcFidHit::Class_Name()
{
   return "tpcFidHit";
}

//______________________________________________________________________________
const char *tpcFidHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::tpcFidHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int tpcFidHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::tpcFidHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *tpcFidHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::tpcFidHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *tpcFidHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::tpcFidHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr scintHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *scintHit::Class_Name()
{
   return "scintHit";
}

//______________________________________________________________________________
const char *scintHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::scintHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int scintHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::scintHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *scintHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::scintHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *scintHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::scintHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr MINDHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MINDHit::Class_Name()
{
   return "MINDHit";
}

//______________________________________________________________________________
const char *MINDHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MINDHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MINDHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MINDHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MINDHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MINDHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MINDHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MINDHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SimulData::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SimulData::Class_Name()
{
   return "SimulData";
}

//______________________________________________________________________________
const char *SimulData::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SimulData*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SimulData::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SimulData*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SimulData::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SimulData*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SimulData::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SimulData*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr GeantBasicParticle::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *GeantBasicParticle::Class_Name()
{
   return "GeantBasicParticle";
}

//______________________________________________________________________________
const char *GeantBasicParticle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GeantBasicParticle*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int GeantBasicParticle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GeantBasicParticle*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GeantBasicParticle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GeantBasicParticle*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GeantBasicParticle::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GeantBasicParticle*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr GeantParticle::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *GeantParticle::Class_Name()
{
   return "GeantParticle";
}

//______________________________________________________________________________
const char *GeantParticle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GeantParticle*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int GeantParticle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GeantParticle*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GeantParticle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GeantParticle*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GeantParticle::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GeantParticle*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr GeantTrackingTruth::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *GeantTrackingTruth::Class_Name()
{
   return "GeantTrackingTruth";
}

//______________________________________________________________________________
const char *GeantTrackingTruth::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GeantTrackingTruth*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int GeantTrackingTruth::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GeantTrackingTruth*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GeantTrackingTruth::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GeantTrackingTruth*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GeantTrackingTruth::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GeantTrackingTruth*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr MemTestRecord::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MemTestRecord::Class_Name()
{
   return "MemTestRecord";
}

//______________________________________________________________________________
const char *MemTestRecord::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MemTestRecord*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MemTestRecord::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MemTestRecord*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MemTestRecord::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MemTestRecord*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MemTestRecord::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MemTestRecord*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void namedRecord::Streamer(TBuffer &R__b)
{
   // Stream an object of class namedRecord.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(namedRecord::Class(),this);
   } else {
      R__b.WriteClassBuffer(namedRecord::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_namedRecord(void *p) {
      return  p ? new(p) ::namedRecord : new ::namedRecord;
   }
   static void *newArray_namedRecord(Long_t nElements, void *p) {
      return p ? new(p) ::namedRecord[nElements] : new ::namedRecord[nElements];
   }
   // Wrapper around operator delete
   static void delete_namedRecord(void *p) {
      delete ((::namedRecord*)p);
   }
   static void deleteArray_namedRecord(void *p) {
      delete [] ((::namedRecord*)p);
   }
   static void destruct_namedRecord(void *p) {
      typedef ::namedRecord current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::namedRecord

//______________________________________________________________________________
void ParticleDescrShortRecord::Streamer(TBuffer &R__b)
{
   // Stream an object of class ParticleDescrShortRecord.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ParticleDescrShortRecord::Class(),this);
   } else {
      R__b.WriteClassBuffer(ParticleDescrShortRecord::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ParticleDescrShortRecord(void *p) {
      return  p ? new(p) ::ParticleDescrShortRecord : new ::ParticleDescrShortRecord;
   }
   static void *newArray_ParticleDescrShortRecord(Long_t nElements, void *p) {
      return p ? new(p) ::ParticleDescrShortRecord[nElements] : new ::ParticleDescrShortRecord[nElements];
   }
   // Wrapper around operator delete
   static void delete_ParticleDescrShortRecord(void *p) {
      delete ((::ParticleDescrShortRecord*)p);
   }
   static void deleteArray_ParticleDescrShortRecord(void *p) {
      delete [] ((::ParticleDescrShortRecord*)p);
   }
   static void destruct_ParticleDescrShortRecord(void *p) {
      typedef ::ParticleDescrShortRecord current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ParticleDescrShortRecord

//______________________________________________________________________________
void ParticleDescrRecord::Streamer(TBuffer &R__b)
{
   // Stream an object of class ParticleDescrRecord.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ParticleDescrRecord::Class(),this);
   } else {
      R__b.WriteClassBuffer(ParticleDescrRecord::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ParticleDescrRecord(void *p) {
      return  p ? new(p) ::ParticleDescrRecord : new ::ParticleDescrRecord;
   }
   static void *newArray_ParticleDescrRecord(Long_t nElements, void *p) {
      return p ? new(p) ::ParticleDescrRecord[nElements] : new ::ParticleDescrRecord[nElements];
   }
   // Wrapper around operator delete
   static void delete_ParticleDescrRecord(void *p) {
      delete ((::ParticleDescrRecord*)p);
   }
   static void deleteArray_ParticleDescrRecord(void *p) {
      delete [] ((::ParticleDescrRecord*)p);
   }
   static void destruct_ParticleDescrRecord(void *p) {
      typedef ::ParticleDescrRecord current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ParticleDescrRecord

//______________________________________________________________________________
void BaseEventRecord::Streamer(TBuffer &R__b)
{
   // Stream an object of class BaseEventRecord.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(BaseEventRecord::Class(),this);
   } else {
      R__b.WriteClassBuffer(BaseEventRecord::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_BaseEventRecord(void *p) {
      return  p ? new(p) ::BaseEventRecord : new ::BaseEventRecord;
   }
   static void *newArray_BaseEventRecord(Long_t nElements, void *p) {
      return p ? new(p) ::BaseEventRecord[nElements] : new ::BaseEventRecord[nElements];
   }
   // Wrapper around operator delete
   static void delete_BaseEventRecord(void *p) {
      delete ((::BaseEventRecord*)p);
   }
   static void deleteArray_BaseEventRecord(void *p) {
      delete [] ((::BaseEventRecord*)p);
   }
   static void destruct_BaseEventRecord(void *p) {
      typedef ::BaseEventRecord current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::BaseEventRecord

//______________________________________________________________________________
void NeutrinoHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class NeutrinoHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(NeutrinoHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(NeutrinoHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_NeutrinoHit(void *p) {
      return  p ? new(p) ::NeutrinoHit : new ::NeutrinoHit;
   }
   static void *newArray_NeutrinoHit(Long_t nElements, void *p) {
      return p ? new(p) ::NeutrinoHit[nElements] : new ::NeutrinoHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_NeutrinoHit(void *p) {
      delete ((::NeutrinoHit*)p);
   }
   static void deleteArray_NeutrinoHit(void *p) {
      delete [] ((::NeutrinoHit*)p);
   }
   static void destruct_NeutrinoHit(void *p) {
      typedef ::NeutrinoHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::NeutrinoHit

//______________________________________________________________________________
void MuonDecayEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class MuonDecayEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MuonDecayEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(MuonDecayEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_MuonDecayEvent(void *p) {
      return  p ? new(p) ::MuonDecayEvent : new ::MuonDecayEvent;
   }
   static void *newArray_MuonDecayEvent(Long_t nElements, void *p) {
      return p ? new(p) ::MuonDecayEvent[nElements] : new ::MuonDecayEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_MuonDecayEvent(void *p) {
      delete ((::MuonDecayEvent*)p);
   }
   static void deleteArray_MuonDecayEvent(void *p) {
      delete [] ((::MuonDecayEvent*)p);
   }
   static void destruct_MuonDecayEvent(void *p) {
      typedef ::MuonDecayEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MuonDecayEvent

//______________________________________________________________________________
void PionDecayEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class PionDecayEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PionDecayEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(PionDecayEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PionDecayEvent(void *p) {
      return  p ? new(p) ::PionDecayEvent : new ::PionDecayEvent;
   }
   static void *newArray_PionDecayEvent(Long_t nElements, void *p) {
      return p ? new(p) ::PionDecayEvent[nElements] : new ::PionDecayEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_PionDecayEvent(void *p) {
      delete ((::PionDecayEvent*)p);
   }
   static void deleteArray_PionDecayEvent(void *p) {
      delete [] ((::PionDecayEvent*)p);
   }
   static void destruct_PionDecayEvent(void *p) {
      typedef ::PionDecayEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PionDecayEvent

//______________________________________________________________________________
void NeutrinoEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class NeutrinoEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(NeutrinoEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(NeutrinoEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_NeutrinoEvent(void *p) {
      return  p ? new(p) ::NeutrinoEvent : new ::NeutrinoEvent;
   }
   static void *newArray_NeutrinoEvent(Long_t nElements, void *p) {
      return p ? new(p) ::NeutrinoEvent[nElements] : new ::NeutrinoEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_NeutrinoEvent(void *p) {
      delete ((::NeutrinoEvent*)p);
   }
   static void deleteArray_NeutrinoEvent(void *p) {
      delete [] ((::NeutrinoEvent*)p);
   }
   static void destruct_NeutrinoEvent(void *p) {
      typedef ::NeutrinoEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::NeutrinoEvent

//______________________________________________________________________________
void SDHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class SDHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SDHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(SDHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SDHit(void *p) {
      return  p ? new(p) ::SDHit : new ::SDHit;
   }
   static void *newArray_SDHit(Long_t nElements, void *p) {
      return p ? new(p) ::SDHit[nElements] : new ::SDHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_SDHit(void *p) {
      delete ((::SDHit*)p);
   }
   static void deleteArray_SDHit(void *p) {
      delete [] ((::SDHit*)p);
   }
   static void destruct_SDHit(void *p) {
      typedef ::SDHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SDHit

//______________________________________________________________________________
void tpcFidHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class tpcFidHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(tpcFidHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(tpcFidHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_tpcFidHit(void *p) {
      return  p ? new(p) ::tpcFidHit : new ::tpcFidHit;
   }
   static void *newArray_tpcFidHit(Long_t nElements, void *p) {
      return p ? new(p) ::tpcFidHit[nElements] : new ::tpcFidHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_tpcFidHit(void *p) {
      delete ((::tpcFidHit*)p);
   }
   static void deleteArray_tpcFidHit(void *p) {
      delete [] ((::tpcFidHit*)p);
   }
   static void destruct_tpcFidHit(void *p) {
      typedef ::tpcFidHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::tpcFidHit

//______________________________________________________________________________
void scintHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class scintHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(scintHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(scintHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_scintHit(void *p) {
      return  p ? new(p) ::scintHit : new ::scintHit;
   }
   static void *newArray_scintHit(Long_t nElements, void *p) {
      return p ? new(p) ::scintHit[nElements] : new ::scintHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_scintHit(void *p) {
      delete ((::scintHit*)p);
   }
   static void deleteArray_scintHit(void *p) {
      delete [] ((::scintHit*)p);
   }
   static void destruct_scintHit(void *p) {
      typedef ::scintHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::scintHit

//______________________________________________________________________________
void MINDHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class MINDHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MINDHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(MINDHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_MINDHit(void *p) {
      return  p ? new(p) ::MINDHit : new ::MINDHit;
   }
   static void *newArray_MINDHit(Long_t nElements, void *p) {
      return p ? new(p) ::MINDHit[nElements] : new ::MINDHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_MINDHit(void *p) {
      delete ((::MINDHit*)p);
   }
   static void deleteArray_MINDHit(void *p) {
      delete [] ((::MINDHit*)p);
   }
   static void destruct_MINDHit(void *p) {
      typedef ::MINDHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MINDHit

//______________________________________________________________________________
void SimulData::Streamer(TBuffer &R__b)
{
   // Stream an object of class SimulData.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SimulData::Class(),this);
   } else {
      R__b.WriteClassBuffer(SimulData::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SimulData(void *p) {
      return  p ? new(p) ::SimulData : new ::SimulData;
   }
   static void *newArray_SimulData(Long_t nElements, void *p) {
      return p ? new(p) ::SimulData[nElements] : new ::SimulData[nElements];
   }
   // Wrapper around operator delete
   static void delete_SimulData(void *p) {
      delete ((::SimulData*)p);
   }
   static void deleteArray_SimulData(void *p) {
      delete [] ((::SimulData*)p);
   }
   static void destruct_SimulData(void *p) {
      typedef ::SimulData current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SimulData

//______________________________________________________________________________
void GeantBasicParticle::Streamer(TBuffer &R__b)
{
   // Stream an object of class GeantBasicParticle.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(GeantBasicParticle::Class(),this);
   } else {
      R__b.WriteClassBuffer(GeantBasicParticle::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_GeantBasicParticle(void *p) {
      return  p ? new(p) ::GeantBasicParticle : new ::GeantBasicParticle;
   }
   static void *newArray_GeantBasicParticle(Long_t nElements, void *p) {
      return p ? new(p) ::GeantBasicParticle[nElements] : new ::GeantBasicParticle[nElements];
   }
   // Wrapper around operator delete
   static void delete_GeantBasicParticle(void *p) {
      delete ((::GeantBasicParticle*)p);
   }
   static void deleteArray_GeantBasicParticle(void *p) {
      delete [] ((::GeantBasicParticle*)p);
   }
   static void destruct_GeantBasicParticle(void *p) {
      typedef ::GeantBasicParticle current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::GeantBasicParticle

//______________________________________________________________________________
void GeantParticle::Streamer(TBuffer &R__b)
{
   // Stream an object of class GeantParticle.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(GeantParticle::Class(),this);
   } else {
      R__b.WriteClassBuffer(GeantParticle::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_GeantParticle(void *p) {
      return  p ? new(p) ::GeantParticle : new ::GeantParticle;
   }
   static void *newArray_GeantParticle(Long_t nElements, void *p) {
      return p ? new(p) ::GeantParticle[nElements] : new ::GeantParticle[nElements];
   }
   // Wrapper around operator delete
   static void delete_GeantParticle(void *p) {
      delete ((::GeantParticle*)p);
   }
   static void deleteArray_GeantParticle(void *p) {
      delete [] ((::GeantParticle*)p);
   }
   static void destruct_GeantParticle(void *p) {
      typedef ::GeantParticle current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::GeantParticle

//______________________________________________________________________________
void GeantTrackingTruth::Streamer(TBuffer &R__b)
{
   // Stream an object of class GeantTrackingTruth.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(GeantTrackingTruth::Class(),this);
   } else {
      R__b.WriteClassBuffer(GeantTrackingTruth::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_GeantTrackingTruth(void *p) {
      return  p ? new(p) ::GeantTrackingTruth : new ::GeantTrackingTruth;
   }
   static void *newArray_GeantTrackingTruth(Long_t nElements, void *p) {
      return p ? new(p) ::GeantTrackingTruth[nElements] : new ::GeantTrackingTruth[nElements];
   }
   // Wrapper around operator delete
   static void delete_GeantTrackingTruth(void *p) {
      delete ((::GeantTrackingTruth*)p);
   }
   static void deleteArray_GeantTrackingTruth(void *p) {
      delete [] ((::GeantTrackingTruth*)p);
   }
   static void destruct_GeantTrackingTruth(void *p) {
      typedef ::GeantTrackingTruth current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::GeantTrackingTruth

//______________________________________________________________________________
void MemTestRecord::Streamer(TBuffer &R__b)
{
   // Stream an object of class MemTestRecord.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MemTestRecord::Class(),this);
   } else {
      R__b.WriteClassBuffer(MemTestRecord::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_MemTestRecord(void *p) {
      return  p ? new(p) ::MemTestRecord : new ::MemTestRecord;
   }
   static void *newArray_MemTestRecord(Long_t nElements, void *p) {
      return p ? new(p) ::MemTestRecord[nElements] : new ::MemTestRecord[nElements];
   }
   // Wrapper around operator delete
   static void delete_MemTestRecord(void *p) {
      delete ((::MemTestRecord*)p);
   }
   static void deleteArray_MemTestRecord(void *p) {
      delete [] ((::MemTestRecord*)p);
   }
   static void destruct_MemTestRecord(void *p) {
      typedef ::MemTestRecord current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MemTestRecord

namespace ROOT {
   static TClass *vectorlEscintHitgR_Dictionary();
   static void vectorlEscintHitgR_TClassManip(TClass*);
   static void *new_vectorlEscintHitgR(void *p = 0);
   static void *newArray_vectorlEscintHitgR(Long_t size, void *p);
   static void delete_vectorlEscintHitgR(void *p);
   static void deleteArray_vectorlEscintHitgR(void *p);
   static void destruct_vectorlEscintHitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<scintHit>*)
   {
      vector<scintHit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<scintHit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<scintHit>", -2, "vector", 214,
                  typeid(vector<scintHit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEscintHitgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<scintHit>) );
      instance.SetNew(&new_vectorlEscintHitgR);
      instance.SetNewArray(&newArray_vectorlEscintHitgR);
      instance.SetDelete(&delete_vectorlEscintHitgR);
      instance.SetDeleteArray(&deleteArray_vectorlEscintHitgR);
      instance.SetDestructor(&destruct_vectorlEscintHitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<scintHit> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<scintHit>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEscintHitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<scintHit>*)0x0)->GetClass();
      vectorlEscintHitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEscintHitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEscintHitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<scintHit> : new vector<scintHit>;
   }
   static void *newArray_vectorlEscintHitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<scintHit>[nElements] : new vector<scintHit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEscintHitgR(void *p) {
      delete ((vector<scintHit>*)p);
   }
   static void deleteArray_vectorlEscintHitgR(void *p) {
      delete [] ((vector<scintHit>*)p);
   }
   static void destruct_vectorlEscintHitgR(void *p) {
      typedef vector<scintHit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<scintHit>

namespace ROOT {
   static TClass *vectorlESDHitgR_Dictionary();
   static void vectorlESDHitgR_TClassManip(TClass*);
   static void *new_vectorlESDHitgR(void *p = 0);
   static void *newArray_vectorlESDHitgR(Long_t size, void *p);
   static void delete_vectorlESDHitgR(void *p);
   static void deleteArray_vectorlESDHitgR(void *p);
   static void destruct_vectorlESDHitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<SDHit>*)
   {
      vector<SDHit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<SDHit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<SDHit>", -2, "vector", 214,
                  typeid(vector<SDHit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlESDHitgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<SDHit>) );
      instance.SetNew(&new_vectorlESDHitgR);
      instance.SetNewArray(&newArray_vectorlESDHitgR);
      instance.SetDelete(&delete_vectorlESDHitgR);
      instance.SetDeleteArray(&deleteArray_vectorlESDHitgR);
      instance.SetDestructor(&destruct_vectorlESDHitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<SDHit> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<SDHit>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlESDHitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<SDHit>*)0x0)->GetClass();
      vectorlESDHitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlESDHitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlESDHitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<SDHit> : new vector<SDHit>;
   }
   static void *newArray_vectorlESDHitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<SDHit>[nElements] : new vector<SDHit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlESDHitgR(void *p) {
      delete ((vector<SDHit>*)p);
   }
   static void deleteArray_vectorlESDHitgR(void *p) {
      delete [] ((vector<SDHit>*)p);
   }
   static void destruct_vectorlESDHitgR(void *p) {
      typedef vector<SDHit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<SDHit>

namespace ROOT {
   static TClass *vectorlEParticleDescrShortRecordgR_Dictionary();
   static void vectorlEParticleDescrShortRecordgR_TClassManip(TClass*);
   static void *new_vectorlEParticleDescrShortRecordgR(void *p = 0);
   static void *newArray_vectorlEParticleDescrShortRecordgR(Long_t size, void *p);
   static void delete_vectorlEParticleDescrShortRecordgR(void *p);
   static void deleteArray_vectorlEParticleDescrShortRecordgR(void *p);
   static void destruct_vectorlEParticleDescrShortRecordgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ParticleDescrShortRecord>*)
   {
      vector<ParticleDescrShortRecord> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ParticleDescrShortRecord>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ParticleDescrShortRecord>", -2, "vector", 214,
                  typeid(vector<ParticleDescrShortRecord>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEParticleDescrShortRecordgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<ParticleDescrShortRecord>) );
      instance.SetNew(&new_vectorlEParticleDescrShortRecordgR);
      instance.SetNewArray(&newArray_vectorlEParticleDescrShortRecordgR);
      instance.SetDelete(&delete_vectorlEParticleDescrShortRecordgR);
      instance.SetDeleteArray(&deleteArray_vectorlEParticleDescrShortRecordgR);
      instance.SetDestructor(&destruct_vectorlEParticleDescrShortRecordgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ParticleDescrShortRecord> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<ParticleDescrShortRecord>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEParticleDescrShortRecordgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<ParticleDescrShortRecord>*)0x0)->GetClass();
      vectorlEParticleDescrShortRecordgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEParticleDescrShortRecordgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEParticleDescrShortRecordgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<ParticleDescrShortRecord> : new vector<ParticleDescrShortRecord>;
   }
   static void *newArray_vectorlEParticleDescrShortRecordgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<ParticleDescrShortRecord>[nElements] : new vector<ParticleDescrShortRecord>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEParticleDescrShortRecordgR(void *p) {
      delete ((vector<ParticleDescrShortRecord>*)p);
   }
   static void deleteArray_vectorlEParticleDescrShortRecordgR(void *p) {
      delete [] ((vector<ParticleDescrShortRecord>*)p);
   }
   static void destruct_vectorlEParticleDescrShortRecordgR(void *p) {
      typedef vector<ParticleDescrShortRecord> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<ParticleDescrShortRecord>

namespace ROOT {
   static TClass *vectorlEGeantParticlegR_Dictionary();
   static void vectorlEGeantParticlegR_TClassManip(TClass*);
   static void *new_vectorlEGeantParticlegR(void *p = 0);
   static void *newArray_vectorlEGeantParticlegR(Long_t size, void *p);
   static void delete_vectorlEGeantParticlegR(void *p);
   static void deleteArray_vectorlEGeantParticlegR(void *p);
   static void destruct_vectorlEGeantParticlegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<GeantParticle>*)
   {
      vector<GeantParticle> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<GeantParticle>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<GeantParticle>", -2, "vector", 214,
                  typeid(vector<GeantParticle>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEGeantParticlegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<GeantParticle>) );
      instance.SetNew(&new_vectorlEGeantParticlegR);
      instance.SetNewArray(&newArray_vectorlEGeantParticlegR);
      instance.SetDelete(&delete_vectorlEGeantParticlegR);
      instance.SetDeleteArray(&deleteArray_vectorlEGeantParticlegR);
      instance.SetDestructor(&destruct_vectorlEGeantParticlegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<GeantParticle> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<GeantParticle>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEGeantParticlegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<GeantParticle>*)0x0)->GetClass();
      vectorlEGeantParticlegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEGeantParticlegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEGeantParticlegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<GeantParticle> : new vector<GeantParticle>;
   }
   static void *newArray_vectorlEGeantParticlegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<GeantParticle>[nElements] : new vector<GeantParticle>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEGeantParticlegR(void *p) {
      delete ((vector<GeantParticle>*)p);
   }
   static void deleteArray_vectorlEGeantParticlegR(void *p) {
      delete [] ((vector<GeantParticle>*)p);
   }
   static void destruct_vectorlEGeantParticlegR(void *p) {
      typedef vector<GeantParticle> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<GeantParticle>

namespace {
  void TriggerDictionaryInitialization_GasTPCDataLibDict_Impl() {
    static const char* headers[] = {
"GasTPCDataLib.hxx",
0
    };
    static const char* includePaths[] = {
"/storage/epp2/epp_software/2015.1/Cellar/root6/6.06.06/include/root",
"/home/epp/phsmaj/hptrex/src/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "GasTPCDataLibDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$GasTPCDataLib.hxx")))  namedRecord;
class __attribute__((annotate("$clingAutoload$GasTPCDataLib.hxx")))  ParticleDescrShortRecord;
class __attribute__((annotate("$clingAutoload$GasTPCDataLib.hxx")))  ParticleDescrRecord;
class __attribute__((annotate("$clingAutoload$GasTPCDataLib.hxx")))  BaseEventRecord;
class __attribute__((annotate("$clingAutoload$GasTPCDataLib.hxx")))  NeutrinoHit;
class __attribute__((annotate("$clingAutoload$GasTPCDataLib.hxx")))  MuonDecayEvent;
class __attribute__((annotate("$clingAutoload$GasTPCDataLib.hxx")))  PionDecayEvent;
class __attribute__((annotate("$clingAutoload$GasTPCDataLib.hxx")))  NeutrinoEvent;
class __attribute__((annotate("$clingAutoload$GasTPCDataLib.hxx")))  SDHit;
class __attribute__((annotate("$clingAutoload$GasTPCDataLib.hxx")))  tpcFidHit;
class __attribute__((annotate("$clingAutoload$GasTPCDataLib.hxx")))  scintHit;
class __attribute__((annotate("$clingAutoload$GasTPCDataLib.hxx")))  MINDHit;
class __attribute__((annotate("$clingAutoload$GasTPCDataLib.hxx")))  SimulData;
class __attribute__((annotate("$clingAutoload$GasTPCDataLib.hxx")))  GeantBasicParticle;
class __attribute__((annotate("$clingAutoload$GasTPCDataLib.hxx")))  GeantParticle;
class __attribute__((annotate("$clingAutoload$GasTPCDataLib.hxx")))  GeantTrackingTruth;
class __attribute__((annotate("$clingAutoload$GasTPCDataLib.hxx")))  MemTestRecord;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "GasTPCDataLibDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "GasTPCDataLib.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"BaseEventRecord", payloadCode, "@",
"GeantBasicParticle", payloadCode, "@",
"GeantParticle", payloadCode, "@",
"GeantTrackingTruth", payloadCode, "@",
"MINDHit", payloadCode, "@",
"MemTestRecord", payloadCode, "@",
"MuonDecayEvent", payloadCode, "@",
"NeutrinoEvent", payloadCode, "@",
"NeutrinoHit", payloadCode, "@",
"ParticleDescrRecord", payloadCode, "@",
"ParticleDescrShortRecord", payloadCode, "@",
"PionDecayEvent", payloadCode, "@",
"SDHit", payloadCode, "@",
"SimulData", payloadCode, "@",
"namedRecord", payloadCode, "@",
"scintHit", payloadCode, "@",
"tpcFidHit", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("GasTPCDataLibDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_GasTPCDataLibDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_GasTPCDataLibDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_GasTPCDataLibDict() {
  TriggerDictionaryInitialization_GasTPCDataLibDict_Impl();
}
