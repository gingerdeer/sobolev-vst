#include "IPlugEffect.h"
#include "IPlug_include_in_plug_src.h"
#include "IControls.h"
#include <cmath>
#include <iostream>
#include <string>
#include <cstdio>
#include "AudioFFT.h"
#include <fft.h>
#include <complex>

IPlugEffect::IPlugEffect(const InstanceInfo& info)
: Plugin(info, MakeConfig(kNumParams, kNumPresets))
{
  GetParam(kGain)->InitDouble("Gain", 100., 0., 300.0, 0.01, "%");
  
  GetParam(kItercount)->InitInt("Itercount", 0, -150, 500.0, "Amount");
  GetParam(kLambda)->InitDouble("Lambda", 1.0, 0.0001, 10.0,0.01, "Awesomeness");
  GetParam(kK)->InitDouble("K", 1.0, 0.0, 100.0, 1.0, "Happy joyness");
  //  bump2 = np.exp(-1.0 * abs(t) / (np.sqrt(llambda)))
  //  bump2 /= np.trapz(bump2) # normalize the integral to 1
  
  //for (int i = 0; i < convLen; ++i) {
   // convArray[i] = exp(-1.0 * abs(start) / sqrt(lambda));
  // / start += step;
 // }
  // approx integral? 0.99?

  




#if IPLUG_EDITOR // http://bit.ly/2S64BDd
  mMakeGraphicsFunc = [&]() {
    return MakeGraphics(*this, PLUG_WIDTH, PLUG_HEIGHT, PLUG_FPS, GetScaleForScreen(PLUG_WIDTH, PLUG_HEIGHT));
  };
  
  mLayoutFunc = [&](IGraphics* pGraphics) {
    pGraphics->AttachCornerResizer(EUIResizerMode::Scale, false);
    const IBitmap bitmap1 = pGraphics->LoadBitmap(BACKGROUND_IM, 1/* num frames*/);
    //pGraphics->AttachPanelBackground(bitmap1);
    pGraphics->AttachBackground(BACKGROUND_IM);
    pGraphics->AttachPanelBackground(COLOR_LIGHT_GRAY);
    pGraphics->LoadFont("Roboto-Regular", ROBOTO_FN);
    const IRECT b = pGraphics->GetBounds();


    pGraphics->AttachControl(new ITextControl(b.GetMidVPadded(50), "Awesomizer v0.1", IText(50, COLOR_INDIGO)),1);
    pGraphics->AttachControl(new IVKnobControl(b.GetCentredInside(100).GetVShifted(100), kGain),2);
    pGraphics->AttachControl(new IVKnobControl(b.GetCentredInside(100).GetVShifted(-200), kItercount),3);
    pGraphics->AttachControl(new IVKnobControl(b.GetCentredInside(100).GetVShifted(-150).GetHShifted(75), kLambda),4);
    pGraphics->AttachControl(new IVKnobControl(b.GetCentredInside(100).GetVShifted(-150).GetHShifted(-75), kK),5);
    pGraphics->AttachControl(new ITextControl(b.GetMidVPadded(50).GetVShifted(200), "Copyright 2021 Pentti Sunila", IText(50, COLOR_INDIGO)));
    pGraphics->AttachControl(new ITextControl(b.GetMidVPadded(50).GetVShifted(250), "All rights reserved", IText(50, COLOR_INDIGO)));
    pGraphics->GetControlWithTag(2)->SetText(IText(14.F, COLOR_BLUE));
    pGraphics->GetControlWithTag(3)->SetText(IText(14.F, COLOR_BLUE));
    pGraphics->GetControlWithTag(4)->SetText(IText(14.F, COLOR_BLUE));
    pGraphics->GetControlWithTag(5)->SetText(IText(14.F, COLOR_BLUE));
    //pGraphics->
  };
#endif
}



#if IPLUG_DSP
void IPlugEffect::ProcessBlock(sample** inputs, sample** outputs, int nFrames)
{
  const double gain = GetParam(kGain)->Value() / 100.;
  const int nChans = NOutChansConnected();
  const int itc = GetParam(kItercount)->Value();


  const double llambda = GetParam(kLambda)->Value();
  const double k = GetParam(kK)->Value();
  audiofft::AudioFFT fft;

  int sz = 2;
  while (sz < nFrames * 2) sz *= 2;
  const size_t fftSize = sz; // Needs to be power of 2!

  fft.init(fftSize);

  // std::
  //double step = 30. / audiofft::AudioFFT::ComplexSize(fftSize);
  std::vector<double> convfilter(audiofft::AudioFFT::ComplexSize(fftSize),1.0f);
  double val = -10.;

  // for (int i = 0; i < convfilter.size(); ++i) {
  //   convfilter[i] = val;
  //   val += step;
  // }



  std::vector<double> convphase(audiofft::AudioFFT::ComplexSize(fftSize), 0.0f);
  double t;
  double step = (2 * abs(val)) / convfilter.size();
  for (int i = 0; i < convfilter.size(); ++i) {
    t = val;
    //convfilter[i] = 1;
    //convfilter[i] = (1.0/sqrt(2.0* M_PI)) *   exp(itc * ((-abs(t) * abs(t)) / pow((1.0 + llambda * abs(t) * abs(t)),  k)      ));
    convfilter[i] =  exp(itc * ((-abs(t) * abs(t)) / pow((1.0 + llambda * abs(t) * abs(t)), k)));
    //convfilter[i] = convfilter[i] < 0.000001 ? 0.0 : convfilter[i];
    //convphase[i] = convfilter[i]; //(1 / sqrt(2 * 3.1415)) * exp(itc * ((-abs(t) * abs(t)) / (1 + abs(t) * abs(t))));
    val += step;
  }

  std::vector<double> inputl(fftSize, 0.0f);
  std::vector<double> inputr(fftSize, 0.0f);

  //for (int i = 0; i < nFrames; ++i) {
  int indx = 0;
  for (int i = 0; i < nFrames * 2; i += 2) {
    inputl[i] = (double)inputs[0][indx % nFrames] * pow(cos(indx * M_PI / nFrames), 2.);
    inputr[i] = (double)inputs[1][indx % nFrames] * pow(cos(indx * M_PI / nFrames), 2.);
   // inputl[i] *= sin(i * M_PI / (nFrames * 2.0)) * sin(i * M_PI / (nFrames * 2.0));
   // inputr[i] *= sin(i * M_PI / (nFrames * 2.0)) * sin(i * M_PI / (nFrames * 2.0));
    indx++;
    inputl[i + 1] = 0; // (double)inputs[0][i % nFrames];
    inputr[i + 1] = 0; // (double)inputs[1][i % nFrames];
  }


  // lmax rmax
  double lmax = *max_element(inputl.begin(),inputl.end());
  double rmax = *max_element(inputr.begin(), inputr.end());


  // normalize?????
  for (int i = 0; i < fftSize; ++i) {
    inputl[i] = inputl[i] / lmax;
    inputr[i] = inputr[i] / rmax;
    //inputl[i] *= sin(i * M_PI / (fftSize)) * sin(i * M_PI / (fftSize));
    //inputr[i] *= sin(i * M_PI / (fftSize)) * sin(i * M_PI / (fftSize));
  }



  std::vector<double> re(audiofft::AudioFFT::ComplexSize(fftSize));
  std::vector<double> im(audiofft::AudioFFT::ComplexSize(fftSize));
  std::vector<double> re2(audiofft::AudioFFT::ComplexSize(fftSize));
  std::vector<double> im2(audiofft::AudioFFT::ComplexSize(fftSize));
  std::vector<double> output(fftSize);
  std::vector<double> output2(fftSize);

  std::vector<double> convput(fftSize);
  //  fft.ifft(convput.data(), convfilter.data(), convphase.data());

  std::vector<double> rec(audiofft::AudioFFT::ComplexSize(fftSize));
  std::vector<double> imc(audiofft::AudioFFT::ComplexSize(fftSize));

  std::vector<double> rec2(audiofft::AudioFFT::ComplexSize(fftSize));
  std::vector<double> imc2(audiofft::AudioFFT::ComplexSize(fftSize));

  fft.init(fftSize);
  fft.fft(inputl.data(), re.data(), im.data());
  fft.fft(inputr.data(), re2.data(), im2.data());

  for (int i = 0; i < re.size(); ++i) {
    rec[i] = re[i] * convfilter[i] - im[i] * convphase[i];
    imc[i] = re[i] * convphase[i] + im[i] * convfilter[i];
    rec2[i] = re2[i] * convfilter[i] - im2[i] * convphase[i];
    imc2[i] = re2[i] * convphase[i] + im2[i] * convfilter[i];
  }

  fft.ifft(output.data(), rec.data(), imc.data());

  fft.ifft(output2.data(), rec2.data(), imc2.data());

  // newlmax newrmax

  // lmax rmax
  double newlmax = *max_element(output.begin(), output.end());
  double newrmax = *max_element(output2.begin(), output2.end());
  // scale by gain* oldmax/newmax
  double scalefactorl = lmax / newlmax;
  double scalefactorr = rmax / newrmax;
  int ind = 0;
  for (int i = 0; i < nFrames * 2; i += 2) {
    outputs[0][ind] = (sample)output.data()[i] * scalefactorl *  gain;
    outputs[1][ind] = (sample)output2.data()[i] * scalefactorr *  gain;
    ind++;
  }

  
  /*

for (int s = 0; s < nFrames; s++) {
  for (int c = 0; c < nChans; c++) {
    outputs[c][s] = (sample)output[c][s] * gain;
  }
}
*/

  /*
  for (int s = 0; s < nFrames; s++) {
    for (int c = 0; c < nChans; c++) {
        for (int i = 0; i < nFrames; ++i) {
          outputs[c][s] += inputs[c][i] * convput[63 - (i % 64)];
        }
        outputs[c][s] = inputs[c][s] * gain;
      }
  }
}
*/
}

#endif
