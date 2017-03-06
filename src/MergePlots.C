#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include "TKey.h"

void MergePlots(){
  TFile fUnmerged("plots_unmerged.root");
  TFile fMerged("plots_merged.root");
  TFile fOut("Displays.root","RECREATE");

  TIter next(fUnmerged.GetListOfKeys());
  TKey *key;

  int i=0;
  while ((key = (TKey*)next())){
    char keyName[20];
    strcpy(keyName,key->GetName());
    if(strstr(keyName,"xy")!=NULL){
      std::cout<<"Drawing for "<<keyName<<std::endl;
      TCanvas * ca=new TCanvas(keyName);
      ca->Divide(2,2);
      TCanvas* cPtr=(TCanvas*)fUnmerged.Get(keyName);
      //cPtr->Draw();
      ca->cd(1);
      cPtr->DrawClonePad();
      cPtr=(TCanvas*)fMerged.Get(keyName);
      ca->cd(3);
      cPtr->DrawClonePad();
      keyName[strlen(keyName)-1]='z';
      cPtr=(TCanvas*)fUnmerged.Get(keyName);
      ca->cd(2);
      cPtr->DrawClonePad();
      cPtr=(TCanvas*)fMerged.Get(keyName);
      ca->cd(4);
      cPtr->DrawClonePad();
      char fileName[20];
      sprintf(fileName,"%s.pdf",keyName);
      fOut.cd();
      ca->Write();
      delete ca;
      if(i++>40)break;
    }
  }
}
