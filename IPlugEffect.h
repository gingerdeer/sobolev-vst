#pragma once

#include "IPlug_include_in_plug_hdr.h"

const int kNumPresets = 1;

enum EParams
{
  kGain = 0,
  kItercount,
  kLambda,
  kK,
  kNumParams
};

using namespace iplug;
using namespace igraphics;

class IPlugEffect final : public Plugin
{
public:
  IPlugEffect(const InstanceInfo& info);
  const float convArray[64] = { 0.00000007861233026453508172429847,
0.00000013071287672349143240191674,
0.00000021734320918659662666171442,
0.00000036138804197122446436438083,
0.00000060089899918460311091213469,
0.00000099914652751517602452893414,
0.00000166133374294228972593715291,
0.00000276238742710009709823456951,
0.00000459316758587366193178087165,
0.00000763730252496403865794180610,
0.00001269894659128304008615262949,
0.00002111520972242986941662390432,
0.00003510937528694266804383389280,
0.00005837821405723310941478232383,
0.00009706854219589553466572551255,
0.00016140099583722913195206472015,
0.00026836996691138924418101874281,
0.00044623292914904868176895202936,
0.00074197507772055164374652713732,
0.00123372118012055791494430145860,
0.00205137341668410865674898602151,
0.00341092700886193572057392131569,
0.00567152862816657984956769666951,
0.00943035042864945154128530901971,
0.01568034211543379782827223323238,
0.02607253364732694395078382854081,
0.04335217980492318129437379070623,
0.07208396081717476089334439848244,
0.11985781177494364557745143429202,
0.19929391893314210570942179856502,
0.19929391893314243877632918611198,
0.11985781177494385374426855150887,
0.07208396081717488579343466881255,
0.04335217980492318129437379070623,
0.02607253364732694395078382854081,
0.01568034211543379782827223323238,
0.00943035042864945154128530901971,
0.00567152862816660066624940839120,
0.00341092700886194179210608723452,
0.00205137341668411212619593797513,
0.00123372118012056008334864642961,
0.00074197507772055164374652713732,
0.00044623292914904868176895202936,
0.00026836996691138924418101874281,
0.00016140099583722970115820527504,
0.00009706854219589579216374147785,
0.00005837821405723321105873599435,
0.00003510937528694272903020609511,
0.00002111520972242990329794179449,
0.00001269894659128304008615262949,
0.00000763730252496403865794180610,
0.00000459316758587366193178087165,
0.00000276238742710009709823456951,
0.00000166133374294229565516778369,
0.00000099914652751517771859482864,
0.00000060089899918460406382420035,
0.00000036138804197122515257865048,
0.00000021734320918659662666171442,
0.00000013071287672349143240191674,
0.00000007861233026453508172429847
  };

#if IPLUG_DSP // http://bit.ly/2S64BDd
  void ProcessBlock(sample** inputs, sample** outputs, int nFrames) override;
#endif
};