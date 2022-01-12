////////////////////////////////////////////////
//Author: Sophia Devinyak (s.devinyak@gmail.com , sdevinya@uwaterloo.ca)
//Date created: November 2021 at TRIUMF under the supervision of Dr. Roger Caballero-Folch



#ifndef __input_h__
#define __input_h__

#include <TH1D.h>
#include <TH1F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <vector>
#include <algorithm>
#include <string.h>
#include <TSpectrum.h>
#include <TObject.h>
#include <functional>
#include <iostream>
#include <fstream>

#include <RQ_OBJECT.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGListBox.h>
#include <TGNumberEntry.h>
#include <TDatime.h>
#include <TApplication.h>

////////////////////////////////////////////////////////////////////////////////
/// Creates a window with buttons and input fields to enter the information about the .root file with the gamma spectum of the calibration run.
class input
{
        RQ_OBJECT("input")
        private:
        TGMainFrame *fMain;   				///< The main frame of the window.
        TGTextEntry *f0;				///< Text entry for the path to the .root file with the gamma spectrum histogram.
        TGTextEntry *f2;				///< Text entry for the canvas name in that .root file that contains the histogram.
        TGTextEntry *f3;				///< Text entry for the histogram name.
        TGListBox *fListBox;				///< Listbox to select the radioactive calibration source used in this histogram.
        TGNumberEntry *fN2;				///< Number entry for the reference activity of the source.
        TGNumberEntry *fN3;				///< Number entry for the date when the reference activity of the source was measured.
        TGNumberEntry *fN4;				///< Number entry for the date of the calibration run.
        TGNumberEntry *fN5;				///< Number entry for the length of the calibration run.
        TGNumberEntry *fN6;				///< Number entry for the uncertainty of the calibration run length.
        TGNumberEntry *fN7;				///< Number entry for the uncertainty in the reference activity of the source.
        public:
        input(const TGWindow *p, UInt_t w, UInt_t h);	///< Class constructor: constructs the window for inputting all the information about the file and the run.
        virtual ~input();				///< Class destructor.
        void searchGamma();				///< Process the inputted information and initiate gamma search.
};
////////////////////////////////////////////////////////////////////////////////
/// This class is for the window with searching for gamma peaks and all related functions.
class gammaSearch
{
        RQ_OBJECT("gammaSearch")
        private:
        TRootEmbeddedCanvas *fEcanvas0;				///< Embedded canvas used to plot the gamm aspectrum.
        TGNumberEntry *fN0;					///< Number entry field to enter the sensitivity of the search function.
        TGNumberEntry *fN1;					///< Number entry field to select a peak.
        TGraph *gscat;						///< Scatter plot that will contain the markers that will mark all the found peaks.
        TGraph *selected;					///< Scatter plot that will contain the marker that will mark the selected peak.
        double lastPar=0.0;					///< The initial value of the sensitivity parameter.
        public:
        TGMainFrame *fMain;					///< The main frame of the window.
        gammaSearch(const TGWindow *p, UInt_t w, UInt_t h);	///< Class constructor: constructs a window for gamma search.
        virtual ~gammaSearch();					///< Class destructor.
        void redoSearch();					///< Deletes the previous gamma search results and does it over again using the new sensitivity value.
        void selectPeak();					///< Selects a gamma peak.
        void deletePeak();					///< Deletes the selected gamma peak.
        void correlatePeaks();					///< Wrapper function for performing the correlation of the found gamma peaks to those from the literature.
        void ratios();						///< This function performs the correlations between the found gamma peaks and those from the literature.
};
////////////////////////////////////////////////////////////////////////////////
/// This class is for correlating the detector channels to energy levels and all related functions.
class ratioPeaks
{
        RQ_OBJECT("ratioPeaks")
        private:
        TGMainFrame *fMain1;					///< The main frame of the window.
        TRootEmbeddedCanvas *fEcanvas1;				///< The embedded window where the gamma spectrum will be plotted again with gamma peaks that passed the ratio test marked.
        TGNumberEntry *fN0;					///< Number entry field to select a peak
        TCanvas *keV;						///< Canvas that will contain the plot of energy (keV) vs detector channels
        TGraph *gscat;						///< Scatter plot that will contain the markers that will mark all the peaks that passed the ratio test.
        TGraph *selected;					///< Scatter plot that will contain the marker that will mark the selected peak.
        public:
        ratioPeaks(const TGWindow *p, UInt_t w, UInt_t h);	///< Class constructor: constructs a window that will display the gamma spectrum and creates the correlation of detector channels to energy levels.
        virtual ~ratioPeaks();					///< Class destructor.
        void selectPeak();					///< Selects a peak.
        void deletePeak();					///< Deletes a selected peak.
        void fitPeaks();					///< Initiates fitting of gaussians to the gamma peaks
        void back();						///< Goes back to the previous (gamma search) class, destructing all the windows created by this (ratioPeaks) class.
};
////////////////////////////////////////////////////////////////////////////////
/// This class is for fitting gaussian curves to the gamma peaks.
class gammaFits
{
        RQ_OBJECT("gammaFits")
        private:
        TGMainFrame *fMain;					///< The main frame of the window.
        TRootEmbeddedCanvas *fEcanvas2;				///< The embedded canvas that will contain the gamma peak fits.
        TGLabel *label;						///< The label that will display the energy of the peak and # of gamma particles detected in that peak.
        TCanvas *gausFit;					///< Alias for fEcanvas2. Needed because this canvas is going to be called outside of the constructor too.
        TH1F *fitting;						///< This is the histogram that will contain the fitted gamma peaks.
	TGNumberEntry *fN0;					///< This number entry is for the user to input a custom value of the # of gamma particles in the peak if the program doesn't fit it properly.
        float area;						///< This value will contain the integral of the gamma peak fit
        float nrg;						///< This value will contain the centroid of the gamma peak fit
        unsigned int peak_index = 0;				///< This value will contain the index of the peak that's being fitted at the moment.
        public:			
        gammaFits(const TGWindow *p, UInt_t w, UInt_t h);	///< Class constructor: constructs a window for gamma fits and fits the first gamma peak.
        virtual ~gammaFits();					///< Class destructor.
        void gyes();						///< Approves the fitted peak.
        void gno();						///< Disaproves the fitted peak; checks if the user has inputted a custom valus for the # of gamma particles in the peak.
        void gnext();						///< Moves onto the next peak. Called by gyes() and gno() functions.
};
////////////////////////////////////////////////////////////////////////////////
/// This class is used to plot and fit the efficiency curve.
class efficiency
{
        RQ_OBJECT("efficiency")
        private:
        TGMainFrame *fMain;					///< The main frame of the window.
        TRootEmbeddedCanvas *fEcanvas3;				///< The embedded canvas that will contain the efficiency curve.
        TCanvas *effi;						///< Alias for fEcanvas3. Needed because this canvas is going to be called outside of the constructor too.
        public:
        efficiency(const TGWindow *p, UInt_t w, UInt_t h);	///< Class constructor: constructs a window for efficiency curve.
        virtual ~efficiency();					///< Class destructor.
        void plot();						///< Plots and fits the efficiency curve.
};
#endif
