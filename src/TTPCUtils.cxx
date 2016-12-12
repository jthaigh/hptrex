#include "TTPCUtils.hxx"
#include "TTRExHVCluster.hxx"
#include "TTPCHitPad.hxx"
#include "TTPCLayout.hxx"

namespace TTPCUtils {

  bool Curvature_to_MomentumAndCharge(const TVector3& pos, const TVector3& dir, double curv, double& p, double& q) {
    //*****************************************************************************

    trex::TTPCLayout layout;

    double B = layout.GetBField();
    // project into the bending plane
    double factor = -(0.3 * B) / sqrt(1. - dir.X() * dir.X());

    if (fabs(curv) > 0 && factor != 0) {
      p = fabs(factor / curv);
      q = -curv / fabs(curv);
      return true;
    }

    return false;
  }

  bool MomentumAndCharge_to_Curvature(const TVector3& pos, const TVector3& dir, double p, double q, double& curv){

  trex::TTPCLayout layout;
  
  double B = layout.GetBField();

  // project into the bending plane
  double factor = -(0.3*B)/sqrt(1.-dir.X()*dir.X());
  if (p > 0){
    curv = q*factor/fabs(p);
    return true;
  }
   
  return false;
}

  //*****************************************************************************
  bool SafeSort( double first, double second ){
    // Use int because I don't trust the ordering of doubles
    // since the bug 871 was discovered.
    // Millimeter precision is enough.
    return ( int(first) < int(second) );
  }

  //*****************************************************************************
  int SenseFromTwoClusters(trex::TTRExPath& Path2, trex::TTRExHVCluster& ChosenClu){
  //*****************************************************************************
    int sense = 0;
    if (Path2.GetClusters().size() < 2)
      return 0;
    trex::TTRExHVCluster* tmpClu;
    trex::TTRExHVCluster* tmprClu;
    trex::TTRExHVCluster* nextClu;
    auto Clu = Path2.GetClusters().begin();
    auto rClu = Path2.GetClusters().rbegin();
    tmpClu = *Clu;
    tmprClu = *rClu;
    if (&ChosenClu == tmpClu){
      Clu++;
      nextClu = *Clu;
    } else if (&ChosenClu == tmprClu){
      rClu++;
      nextClu = *rClu;
    } else {
      // PROBLEM. That will probably be a bigger problem outside of this function regardless of the sense.
      return 0;
    }

    // Find the crude track "sense" in z(y) when the first two clusters are vertical(horizontal).
    // If the first two clusters have different orientation, just use 0 and deal with it later.
    if ( ChosenClu.IsVertical() == nextClu->IsVertical()){
      if (ChosenClu.IsVertical()){
        sense = nextClu->Z() - ChosenClu.Z();
      }else{
        sense = nextClu->Y() - ChosenClu.Y();
      }
    }
    return sense;
  }


  //*****************************************************************************
  void FindClosestEnds(trex::TTRExPath& PathA, trex::TTRExPath& PathB, unsigned int &EndA, unsigned int &EndB){
    // TODO: Add check on the presence of clusters and use exception if not there.
    trex::TTRExHVCluster& FirstCluA = **(PathA.GetClusters().begin());
    trex::TTRExHVCluster& LastCluA = **(PathA.GetClusters().rbegin());
    trex::TTRExHVCluster& FirstCluB = **(PathB.GetClusters().begin());
    trex::TTRExHVCluster& LastCluB = **(PathB.GetClusters().rbegin());
    std::vector<double> Sorter;
    double FtFt = (FirstCluA.GetPosition() - FirstCluB.GetPosition()).Mag();
    double FtLt = (FirstCluA.GetPosition() - LastCluB.GetPosition()).Mag();
    double LtLt = (LastCluA.GetPosition() - LastCluB.GetPosition()).Mag();
    double LtFt = (LastCluA.GetPosition() - FirstCluB.GetPosition()).Mag();
    Sorter.push_back( FtFt );
    Sorter.push_back( FtLt );
    Sorter.push_back( LtLt );
    Sorter.push_back( LtFt );
    std::stable_sort(Sorter.begin(), Sorter.end(), TTPCUtils::SafeSort);
    double MinDist = *(Sorter.begin());
    if (MinDist == FtFt || MinDist == FtLt){
      EndA = 0;
    } else {
      EndA = 1;
    }
    if (MinDist == FtFt || MinDist == LtFt){
      EndB = 0;
    } else {
      EndB = 1;
    }
  }

  //*****************************************************************************
  trex::TTRExPath* MergePaths(trex::TTRExPath& PathA, trex::TTRExPath& PathB){

    //MDH TODO: Warning - we are creating a path here.
    trex::TTRExPath* newPath = new trex::TTRExPath();

    // Create a new path putting in the clusters from the two matched paths.
    // BE CAREFUL with the ordering of the clusters !
    unsigned int UseEndA, UseEndB;
    FindClosestEnds(PathA, PathB, UseEndA, UseEndB);

    std::vector<trex::TTRExHVCluster*>& newClusters = newPath->GetClusters();
    if (UseEndA == 1){
       std::vector<trex::TTRExHVCluster*>& clusters = PathA.GetClusters();
      newClusters = clusters;
      if (UseEndB == 0 ){
        for (auto hit = PathB.GetClusters().begin(); hit != PathB.GetClusters().end(); ++hit)
          newClusters.push_back(*hit);
      } else {
        for (auto hit = PathB.GetClusters().rbegin(); hit != PathB.GetClusters().rend(); ++hit)
          newClusters.push_back(*hit);
      }
    } else {
      if (UseEndB == 1 ){
	std::vector<trex::TTRExHVCluster*>& clusters = PathB.GetClusters();
        newClusters = clusters;
        for (auto hit = PathA.GetClusters().begin(); hit != PathA.GetClusters().end(); ++hit)
          newClusters.push_back(*hit);
      } else {
        for (auto hit = PathB.GetClusters().rbegin(); hit != PathB.GetClusters().rend(); ++hit)
          newClusters.push_back(*hit);
        for (auto hit = PathA.GetClusters().begin(); hit != PathA.GetClusters().end(); ++hit)
          newClusters.push_back(*hit);
      }
    }
    //    newPath->SetClusters(newClusters);

    //MDH TODO: Manage IDs properly (do we use them for anything?)
    //newPath->SetId(ND::tpcCalibration().GetPathId());

    return newPath;
  }

  void ReverseStateSenseAndCharge(std::vector<double>& propagState){

    //Velocity components
    propagState[3]*=-1;
    propagState[4]*=-1;
    propagState[5]*=-1;

    //Charge
    propagState[6]*=-1;
  }

}
