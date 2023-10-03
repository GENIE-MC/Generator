/*
Checks the GetParallel and GetOrthogonal functions in MathUtils.h
*/

#include "Framework/Numerical/MathUtils.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace genie;

void test_parallel_orthogonal_genie() {
    // Define a unit vector along z-axis
    TVector3 unitVec(0., 0.,1.);
    std::cout << "unitVec: " << unitVec.X() << " , " << unitVec.Y() << " , "
    << unitVec.Z() << std::endl;

    std::cout << "***************************************************************"<<std::endl;

    // Test case 1: LorentzVector is parallel to the unit vector
    TLorentzVector lv1(0., 0., 5., 5.);
    std::cout << "LorentzVector - Test 1: " << lv1.X() << " , " << lv1.Y()
    << " , " << lv1.Z() << " , " << lv1.T() << std::endl;
    auto parallel1 = utils::math::GetParallel(lv1, unitVec);
    auto orthogonal1 = utils::math::GetOrthogonal(lv1, unitVec);
    std::cout << "Test 1 - Parallel: " << parallel1.X() << " , "
    << parallel1.Y() << " , " << parallel1.Z() << " , " << parallel1.T()
    << ", Orthogonal: " 
    << orthogonal1.X() << " , " << orthogonal1.Y()
    << " , " << orthogonal1.Z() << " , " << orthogonal1.T() << std::endl;
    std::cout << "***************************************************************"<<std::endl;


    // Test case 2: LorentzVector is orthogonal to the unit vector
    TLorentzVector lv2(3., 4., 0., 5.);
    std::cout << "LorentzVector - Test 2: " << lv2.X() << " , " << lv2.Y()
    << " , " << lv2.Z() << " , " << lv2.T() << std::endl;
    auto parallel2 = utils::math::GetParallel(lv2, unitVec);
    auto orthogonal2 = utils::math::GetOrthogonal(lv2, unitVec);
    std::cout << "Test 2 - Parallel: " << parallel2.X() << " , "
    << parallel2.Y() << " , " << parallel2.Z() << " , " << parallel2.T()
    << ", Orthogonal: " 
    << orthogonal2.X() << " , " << orthogonal2.Y()
    << " , " << orthogonal2.Z() << " , " << orthogonal2.T() << std::endl;
    std::cout << "***************************************************************"<<std::endl;


    // Test case 3: LorentzVector is a mix of parallel and orthogonal to the unit vector
    TLorentzVector lv3(3., 4., 5., 7.07107);
    std::cout << "LorentzVector - Test 3: " << lv3.X() << " , " << lv3.Y()
    << " , " << lv3.Z() << " , " << lv3.T() << std::endl;
    auto parallel3 = utils::math::GetParallel(lv3, unitVec);
    auto orthogonal3 = utils::math::GetOrthogonal(lv3, unitVec);
    std::cout << "Test 3 - Parallel: " << parallel3.X() << " , "
    << parallel3.Y() << " , " << parallel3.Z() << " , " << parallel3.T()
    << ", Orthogonal: " 
    << orthogonal3.X() << " , " << orthogonal3.Y()
    << " , " << orthogonal3.Z() << " , " << orthogonal3.T() << std::endl;
}
