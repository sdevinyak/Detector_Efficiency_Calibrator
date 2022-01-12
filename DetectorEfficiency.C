////////////////////////////////////////////////
//Author: Sophia Devinyak (s.devinyak@gmail.com , sdevinya@uwaterloo.ca)
//Date created: November 2021 at TRIUMF under the supervision of Dr. Roger Caballero-Folch


#include "input.h"
#include <TLatex.h>
#include <TFile.h>
#include <TF1.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TGraphErrors.h>

using namespace std;

void ScaleAxis(TAxis *a, std::function<Double_t(Double_t)> Scale);   ///< Used to scale the X axis 
Double_t calibration(Double_t *x, Double_t *par);                    ///< This is a linear function used to correlate the detector channels to energy levels
Double_t gausbkg(Double_t *x, Double_t *par);                        ///< Gauss function used to fit gamma peaks together with bkg 
Double_t bkg(Double_t *x, Double_t *par);                            ///< Used by the gausbkg to fit background when fitting gamma peaks. Uses an error function for better estimate of the background "step" behaviour
Double_t effFunc(Double_t *x, Double_t *par);                        ///< This is the efficiency function used to fit the efficiency curve

//for file retrieving
TCanvas *canv;                                                       ///< This variable is used to store the canvas from the user's .root file with the gamma spectrum to be analysed 
TH1F *hist;							     ///< This variable is used to store the histogram from the user's .root file with the gamma spectrum to be analysed
TH1F *safeCopy;                                                      ///< This variable is used to clone the hist variable for future use 
TH1F *corrected;						     ///< This variable is used to store the histogram where the x-axis was scaled to represent the energy in keV
float Time;							     ///< This variable is used to store the length of the calibration run in seconds
float dtime;							     ///< This variable is used to store the uncertainty in the calibration run length in seconds
double halfLife;              					     ///< This variable is used to store the halflife of the calibration source in seconds
double activity;                                                     ///< This variable is used to store the activity of the calibration source at the time of the calibration run
float dact;                                                          ///< This variable is used to store the uncertainty in the activity of the calibration source at the time of the calibration run
float drefact = 200; // Bq                                           ///< This is the default uncertainty in the reference activity of the source

TCanvas *gSearch;  						     ///< Canvas for the gamma peaks search
std::vector<TLatex *> text;                                          ///< This is used to label the gamma peaks found with their detector channel and energy level
TCanvas *test; // for E vs Ch plot				     ///< This canvas will contain the plot of Energy vs Channel
int iso;                           				     ///< The isotope used in the calibration

// for gamma peaks search
TSpectrum *spectrum;                                                 
Int_t gfound;    						     ///< Number of found gamma peaks
Double_t *gpeaks;                                                    ///< x-positions of all the peaks (vector)
std::vector<float> gHeight;                                          ///< Height of the found gamma peaks 
std::vector<float> gMean;					     ///< Centroids of the found peaks

//for gamma peak correlation
std::vector<int> found;                                              ///< Found peaks whose ratios matched the rations of the peaks from the literature
std::vector<int> lit;   					     ///< Peaks from the literature that the peaks from the calibration matched with
std::vector<int> height; 					     ///< Height of the gamma peaks that passed the ratio test
double m; 							     ///< Slope of the linear function used to correlate detector channels to energy // y=mx+b 
double b; 							     ///< Offset of the linear function used to correlate detector channels to energy

//for efficiency correlation
std::vector<double> energy;               			     ///< Temporarily stores the energy of the peaks that were identified and approved by the user (centroids of the gaussian fits)
std::vector<double> Energy;					     ///< Energy of the identified and approved peaks from literature that are within 2.3 keV of found gamma peaks
std::vector<double> Area; 					     ///< Areas of the gaussian fits of the gamma peaks approved by the user
std::vector<double> dArea;                                           ///< Uncertainty of the area of the gaussian fit
std::vector<double> dEnergy;					     ///< Uncertainty in the energy of the identified and approved peak
std::vector<double> yield;                                           ///< Yield from the literature for every identified and approved peak
std::vector<double> dYield;                                          ///< Uncertainty in the yields
std::vector<double> eff;                                             ///< Efficiency of the detector at every gamma peak
std::vector<double> dEff; 					     ///< Uncertainty of the detector at every gamma peak

// energy peaks from literature
std::vector<string> allIsotopes{
"152Eu","60Co"};  //RADIOACTIVE CALIBRATION SOURCES

std::vector<double> allHalfLives {   // in seconds from NNDC   // HALF LIVES
426272112.,166344192.};

std::vector<double> dallHalfLives {  // in seconds     //UNCERTAINTY IN HALF LIVES
283824.,12096.};

std::vector<float> limits{          
0.7,0.5};  // MINIMUM PEAKS THAT NEED TO BE DETECTED; THIS NUMBER IS A PERCENTAGE

std::vector<std::vector <float>> allEnergy{   // GAMMA ENERGY LEVELS
{121.782,344.278,411.116,443.965,778.904,867.373,964.079,1112.069,1212.948,1299.140,1408.005}, // 244.697
{1173.228,1332.490,1460.821,2505.72,2614.532}}; //keV //last 3 peaks in 60Co are natural 40K, sum of 1173 and 1332, and 208Tl

std::vector<std::vector <float>> alldEnergy{   //UNCERTAINTY IN GAMMA ENERGY LEVELS
{0.001,0.001,0.001,0.003,0.002,0.003,0.018,0.003,0.011,0.009,0.003}, //0.001
{0.003,0.004,0.006,0.005,0.013}};//keV

std::vector<std::vector <float>> allYield{    // YIELD OF EACH PEAK LISTED IN allEnergy, IN THE SAME ORDER
{0.286678,0.26558,0.022372,0.031576,0.129603,0.042584,0.146494,0.136855,0.014263,0.016254,0.210692}, //0.076066
{0.998500,0.999826,0.,0.,0.}};

std::vector<std::vector <float>> alldYield{   // UNCERTAINTY IN YIELD
{0.001456,0.005129,0.000246,0.000297,0.001414,0.000274,0.000719,0.000676,0.000093,0.000193,0.001016}, //0.000403
{0.0003,0.000006,0.,0.,0.}};


////////////////////////////////////////////////////////////////////////////////
/// Creates a window with buttons and input fields. This window is designed for a user to input the information about their .root file that contains the gamma spectrum to be analyzed
input::input(const TGWindow *p, UInt_t w, UInt_t h)
{
        // create a main frame
        fMain = new TGMainFrame (p,w,h);

        // create horizontal frames with text entry fields, labels and buttons
        TGHorizontalFrame *hframe = new TGHorizontalFrame(fMain,200,40);   // create a frame
        TGLabel *fL = new TGLabel(hframe,"Enter the path to the .root file containing the histogram");   //create a label
        hframe->AddFrame(fL, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));   // add the label to the window
        f0 = new TGTextEntry(hframe);		// create a text enty
        hframe->AddFrame(f0, new TGLayoutHints(kLHintsCenterX, 5,5,3,4)); // add the text entry
        fMain->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,2,2,2,2)); // add hframe to parent widget. The order in which the frames are added is the order they will appear in the window

        TGHorizontalFrame *hframe2 = new TGHorizontalFrame(fMain,200,40);
        TGLabel *fL2 = new TGLabel(hframe2,"Enter the name of the canvas containing the histogram");
        hframe2->AddFrame(fL2, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
        f2 = new TGTextEntry(hframe2);
        hframe2->AddFrame(f2, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));
        fMain->AddFrame(hframe2, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

        TGHorizontalFrame *hframe3 = new TGHorizontalFrame(fMain,200,40);
        TGLabel *fL3 = new TGLabel(hframe3,"Enter the name of the histogram");
        hframe3->AddFrame(fL3, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
        f3 = new TGTextEntry(hframe3);
        hframe3->AddFrame(f3, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));
        fMain->AddFrame(hframe3, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

	TGHorizontalFrame *hframe5 = new TGHorizontalFrame(fMain,200,40);
	TGLabel *fL4 = new TGLabel(hframe5,"Select the isotope");
	hframe5->AddFrame(fL4, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
	fListBox = new TGListBox(hframe5,90);    // making a listbox of isotopes
	for (unsigned int i=0;i<allIsotopes.size();i++)   // making listbox entries
	{
		fListBox->AddEntry(allIsotopes[i].c_str(),i);
	}
	fListBox->Resize(150,80);
	hframe5->AddFrame(fListBox, new TGLayoutHints(kLHintsTop|kLHintsLeft, 5,5,5,5));
        fMain->AddFrame(hframe5, new TGLayoutHints(kLHintsCenterX,2,2,2,2));


        TGHorizontalFrame *hframe6 = new TGHorizontalFrame(fMain,200,40);
        TGLabel *fL5 = new TGLabel(hframe6,"Reference activity of the source in Bq:");
        hframe6->AddFrame(fL5, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
        fN2 = new TGNumberEntry(hframe6,0,9,100, TGNumberFormat::kNESInteger,     // making an integer number entry field
						 TGNumberFormat::kNEANonNegative,
						 TGNumberFormat::kNELLimitMinMax,
						 0.0000,9999999);
        hframe6->AddFrame(fN2, new TGLayoutHints(kLHintsCenterX, 5,5,5,5));
	TGLabel *fL10 = new TGLabel(hframe6,"+/-"); 
	hframe6->AddFrame(fL10, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
        fN7 = new TGNumberEntry(hframe6,0,9,100, TGNumberFormat::kNESRealFour,     // making a real number entry field with four digits
                                                 TGNumberFormat::kNEANonNegative,
                                                 TGNumberFormat::kNELLimitMinMax,
                                                 0.0000,1000000.0000);
        hframe6->AddFrame(fN7, new TGLayoutHints(kLHintsCenterX, 5,5,5,5));

	
        fMain->AddFrame(hframe6, new TGLayoutHints(kLHintsCenterX,2,2,2,2)); // add hframe to parent widget. The order in which the frames are added is the order they will appear in the window

	TGHorizontalFrame *hframe7 = new TGHorizontalFrame(fMain,200,40);
	fN3 = new TGNumberEntry(hframe7, 1,9,999,TGNumberFormat::kNESDayMYear);//,
                                       //         1,99999);
        TGLabel *fL6 = new TGLabel(hframe7,"Date of the reference activity measurement, D/M/Y");
        hframe7->AddFrame(fL6, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
        hframe7->AddFrame(fN3, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
	fMain->AddFrame(hframe7, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

	TGHorizontalFrame *hframe8 = new TGHorizontalFrame(fMain,200,40);
	TGLabel *fL7 = new TGLabel(hframe8,"Date of the experiment measurement, D/M/Y");
	hframe8->AddFrame(fL7, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
	fN4 = new TGNumberEntry(hframe8, 1,9,999,TGNumberFormat::kNESDayMYear);//,
	hframe8->AddFrame(fN4, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
	fMain->AddFrame(hframe8, new TGLayoutHints(kLHintsCenterX,2,2,2,2)); 

        TGHorizontalFrame *hframe9 = new TGHorizontalFrame(fMain,200,40);
	TGLabel *fL8 = new TGLabel(hframe9,"Length of the experiment measurement in seconds (by default, the uncertainty is 1s)");
	hframe9->AddFrame(fL8, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
	fN5= new TGNumberEntry(hframe9,0,9,100, TGNumberFormat::kNESRealFour,
                                                 TGNumberFormat::kNEANonNegative,
                                                 TGNumberFormat::kNELLimitMinMax,
                                                 0,1000000);
	hframe9->AddFrame(fN5, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
	TGLabel *fL9 = new TGLabel(hframe9,"+/-");
	hframe9->AddFrame(fL9, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
	fN6 = new TGNumberEntry(hframe9,0,9,100, TGNumberFormat::kNESRealFour,
                                                 TGNumberFormat::kNEANonNegative,
                                                 TGNumberFormat::kNELLimitMinMax,
                                                 0,1000000);
        hframe9->AddFrame(fN6, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
        fMain->AddFrame(hframe9, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

        TGHorizontalFrame *hframe4 = new TGHorizontalFrame(fMain,200,40);
        TGTextButton *zSear = new TGTextButton(hframe4, "&Search for Gamma Peaks");   // making a  button
        zSear->Connect("Clicked()","input",this,"searchGamma()");			// connecting that button to a function
        hframe4->AddFrame(zSear, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));		// adding the button to the window
        TGTextButton *exit = new TGTextButton(hframe4, "&Exit", "gApplication->Terminate(0)");
        hframe4->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));
        fMain->AddFrame(hframe4, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
        f0->SetText("60Co_uncalibrated.root");
        f2->SetText("c1");
        f3->SetText("hE");

        fMain->SetCleanup(kDeepCleanup);

        // Set a name to the main frame
        fMain->SetWindowName("Detector Efficiency - File Input");
        // next 3 lines make the widgets visible     
        // Map all subwindows of main frame
        fMain->MapSubwindows();
        // Initialize the layout algorithm
        fMain->Resize(fMain->GetDefaultSize());  // execute all layout specs for widgets before the top-level window is shown on screen
        // map main frame
        fMain->MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// This function takes all the input values in the Detector Efficiency - File Input window and uses them to retrieve the histogram with the gamma spectrum to be analyzed and calculates the activity of the calibration source at the time of the calibration run by using the reference activity and the dates when the reference activity was measured and when the calibration run happened
/// This function also initiates the function that will search for gamma peaks

void input::searchGamma()
{
	// retrieving the info about the file
        const char *totFile = f0->GetText();    
cout<<totFile<<endl;
        const char *canvName = f2->GetText();
cout<<canvName<<endl;
        const char *histName = f3->GetText();
cout<<histName<<endl;
	// retrieving the histogram
        TFile *_file0 = TFile::Open(totFile);
        canv = (TCanvas*)_file0->Get(canvName);
        hist = (TH1F*)canv->GetPrimitive(histName);
	corrected = (TH1F*)hist->Clone();
	corrected->SetDirectory(0);
        safeCopy = (TH1F*)hist->Clone();
	safeCopy->SetDirectory(0);
        _file0->TFile::Close();
	iso = fListBox->GetSelected(); // getting the index of the element in the allIsotopes

	gHeight.clear();    //cleaning all the vectors if the program is used twice without closing
	gMean.clear();
	found.clear();
	lit.clear();
	height.clear();
	energy.clear();
	Energy.clear();
	Area.clear();
	dArea.clear();
	dEnergy.clear();
	yield.clear();
	dYield.clear();
	eff.clear();
	dEff.clear();

	cout<<iso<<endl;
	// preparing all numbers for calculations
	int yE;
	int mE;
	int dE;
	fN4->GetDate(yE,mE,dE);  // date of the experiment measurement
	TDatime *measurement = new TDatime(yE,mE,dE,0,0,0);
cout<<"Measurement YMD"<<yE<<mE<<dE<<endl;
	int yR;
	int mR;
	int dR;
	fN3->GetDate(yR,mR,dR);  // date of the activity reference 
cout<<"Ref YMD "<<yR<<mR<<dR<<endl;
	TDatime *reference = new TDatime(yR,mR,dR,0,0,0);
	UInt_t timeDifference = measurement->Convert() -reference->Convert();
cout<<"Difference: "<<timeDifference<<endl;
	float refactivity = fN2->GetNumberEntry()->GetNumber();
	Time = fN5->GetNumberEntry()->GetNumber();
	float dtimeDiff=86400; // 24 hours in seconds
	float temp = fN6->GetNumberEntry()->GetNumber();
	if (temp==0) {dtime=1; } // 1s 
	else {dtime = temp;}
	halfLife = allHalfLives[iso];
	double for_exponent = timeDifference*TMath::Log(2)/halfLife; 
cout<<"Exp "<<for_exponent<<endl;
	// calculating activity at the time of the calibration run
	activity = refactivity * TMath::Exp(-for_exponent);
cout<<"Activity at the time of the measurement: "<<activity<<endl;
	float dHL = dallHalfLives[iso];
	// calculating uncertainty in activity
	float drefactivity =fN7->GetNumberEntry()->GetNumber();
	// uncertainty in lambda:
	float dLambda = TMath::Log(2)*dHL/halfLife/halfLife;
	// uncertainty in the exponent of labda*time
	float dLT = timeDifference*TMath::Log(2)/halfLife*sqrt(pow(dLambda/(timeDifference*TMath::Log(2)/halfLife),2)+pow(dtimeDiff/timeDifference,2));
	// uncertainty in the exponent
	float dExp = dLT*TMath::Exp(timeDifference*TMath::Log(2)/halfLife*(-1) );
	// uncertainty in the final calculated activity
	dact = activity*sqrt(pow(drefactivity/refactivity,2)+pow(dExp/TMath::Exp((-1)*timeDifference*TMath::Log(2)/halfLife ),2)); 

	new gammaSearch(gClient->GetRoot(),200,200);  // initiate gamma peaks search
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a window with the histogram. Performs the search for gamma peaks and labels them with green triangular markers.

gammaSearch::gammaSearch(const TGWindow *p, UInt_t w, UInt_t h)
{
        // create a main frame
        fMain = new TGMainFrame (p,w,h);
	
	// create a canvas widget
	fEcanvas0 = new TRootEmbeddedCanvas("Ecanvas0",fMain, 800, 600);
	fMain->AddFrame(fEcanvas0, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10,1));

        // create horizontal frames with text entry fields, labels and buttons
        TGHorizontalFrame *hframe = new TGHorizontalFrame(fMain,200,40);
        TGLabel *fL0 = new TGLabel(hframe,"Enter the sensitivity for peak search (smaller number = higher sensitivity)");
        hframe->AddFrame(fL0, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
        fN0 = new TGNumberEntry(hframe,0,9,100, TGNumberFormat::kNESRealFour,
						 TGNumberFormat::kNEANonNegative,
						 TGNumberFormat::kNELLimitMinMax,
						 0.0000,0.99999);
        hframe->AddFrame(fN0, new TGLayoutHints(kLHintsCenterX, 5,5,5,5));
	fN0->Connect("ValueSet(Long_t)","gammaSearch",this,"redoSearch()");
        fMain->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,2,2,2,2)); // add hframe to parent widget. The order in which the frames are added is the order they will appear in the window

	TGHorizontalFrame *hframe2 = new TGHorizontalFrame(fMain,200,40);
	fN1 = new TGNumberEntry(hframe2, 1,9,999,TGNumberFormat::kNESInteger,
                                                TGNumberFormat::kNEANonNegative,
                                                TGNumberFormat::kNELLimitMinMax,
                                                1,99999);
        fN1->Connect("ValueSet(Long_t)","gammaSearch",this,"selectPeak()");
        fN1->GetNumberEntry()->Connect("ReturnPressed()","gammaSearch",this,"selectPeak()");
        TGLabel *fL1 = new TGLabel(hframe2,"Select the number of the peak to be removed (counted along the x-axis)");
        hframe2->AddFrame(fL1, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
        hframe2->AddFrame(fN1, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
	TGTextButton *del = new TGTextButton(hframe2,"&Delete This Peak");
        del->Connect("Clicked()","gammaSearch",this,"deletePeak()");
        hframe2->AddFrame(del, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));

	fMain->AddFrame(hframe2, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

        TGHorizontalFrame *hframe4 = new TGHorizontalFrame(fMain,200,40);

        TGTextButton *gRedo = new TGTextButton(hframe4, "&Correlate Found Peaks with Energies from Literature");
        gRedo->Connect("Clicked()","gammaSearch",this,"correlatePeaks()");
        hframe4->AddFrame(gRedo, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));

        TGTextButton *exit = new TGTextButton(hframe4, "&Exit", "gApplication->Terminate(0)");
        hframe4->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));
        fMain->AddFrame(hframe4, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

        fMain->SetCleanup(kDeepCleanup);

        // Set a name to the main frame
        fMain->SetWindowName("Gamma Search");
        // next 3 lines make the widgets visible     
        // Map all subwindows of main frame
        fMain->MapSubwindows();
        // Initialize the layout algorithm
        fMain->Resize(fMain->GetDefaultSize());  // execute all layout specs for widgets before the top-level window is shown on screen
        // map main frame
        fMain->MapWindow();
	
	lastPar=0.0005; //last par = smallest peak size / biggest peak size
	fN0->SetNumber(lastPar);

	// search for gamma peaks
	gSearch = fEcanvas0->GetCanvas();
	gSearch->SetLogy(1);
	gSearch->SetTitle("Gamma Peaks Search");
	hist->Draw();
	hist->SetTitle("Gamma Peaks Search");
	spectrum = new TSpectrum(50); //create TSpectrum for peak search
	int binmax = hist->GetMaximumBin();  
	double x = hist->GetBinContent(binmax);
	cout <<"Tallest peak " <<x << endl;
	gfound = spectrum->Search(hist,0.9,"",lastPar);   // CHANGE THE LAST ARG IF NOT ALL PEAKS ARE PLOTTED
	gpeaks = spectrum->GetPositionX();   // getting x positions of all peaks
	// sorting the peaks in the order of appearance along the x axis
	int ind[gfound];
	TMath::Sort(gfound,gpeaks,ind,false);
	for(int p=0;p<gfound;p++)  
	{
		float xp = gpeaks[ind[p]];
		Int_t bin = hist->GetXaxis()->FindBin(xp);
		float yp = hist->GetBinContent(bin);	// getting the y value of each peak
		gHeight.push_back(yp);
		gMean.push_back(xp);
	}
	//plotting a scatter plot of all peaks
	gscat = new TGraph(gMean.size(),&gMean[0],&gHeight[0]);
	gscat->SetMarkerStyle(23);
	gscat->SetMarkerColor(3);
	gscat->Draw("p");
	//creating a "selected" peak
	float x1[1]={0};
	x1[0]=gMean[0];
	float y1[1]={0};
	y1[0]=gHeight[0];
	selected = new TGraph(1,x1,y1);
	selected->SetMarkerStyle(23);
	selected->SetMarkerColor(1);
	selected->Draw("p");
	gSearch->Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Deletes the previous gamma search results and does it over again using the new sensitivity value. 
/// This is called when the sensitivity of hte search is changed by using arrows or pressing ENTER after entering the sensitivity
void gammaSearch::redoSearch()
{
	lastPar = fN0->GetNumberEntry()->GetNumber();
	gscat->Delete();
	selected->Delete();
	gHeight.clear();
	gMean.clear();
	spectrum = new TSpectrum(50); //create TSpectrum for peak search
	gfound = spectrum->Search(hist,0.9,"",lastPar);   // CHANGE THE LAST ARG IF NOT ALL PEAKS ARE PLOTTED
	gpeaks = spectrum->GetPositionX();   // getting x positions of all peaks
	// sorting the peaks in the order of appearance along the x axis
	int ind[gfound];
	TMath::Sort(gfound,gpeaks,ind,false);
	for(int p=0;p<gfound;p++)  
	{
		float xp = gpeaks[ind[p]];
		Int_t bin = hist->GetXaxis()->FindBin(xp);
		float yp = hist->GetBinContent(bin);	// getting the y value of each peak
		gHeight.push_back(yp);
		gMean.push_back(xp);
	}
	//plotting a scatter plot of all peaks
	gscat = new TGraph(gMean.size(),&gMean[0],&gHeight[0]);
	gscat->SetMarkerStyle(23);
	gscat->SetMarkerColor(3);
	gscat->Draw("p");
//	for (int p=0;p<gfound;p++) {cout<<gMean[p]<<endl;}
	//creating a "selected" peak
	float x1[1]={0};
	x1[0]=gMean[0];
	float y1[1]={0};
	y1[0]=gHeight[0];
	gSearch->cd();
	selected = new TGraph(1,x1,y1);
	selected->SetMarkerStyle(23);
	selected->SetMarkerColor(1);
	selected->Draw("p");
	gSearch->Update();
	fN1->SetNumber(1);
}
unsigned int gpk =0; ///< This value stores hte index of a selected peak

////////////////////////////////////////////////////////////////////////////////
/// Selects a gamma peak and highlights it with a black triangle marker
/// This is done by using the second number entry. The function is called when the arrows on the number entry are used or when ENTER is pressed after entering the integer value
void gammaSearch::selectPeak()
{
        // get the number entry and convert it to 0 = 1st peak base
        TObject *sender = (TObject *)gTQSender;
        if (sender && sender->InheritsFrom("TGNumberEntry"))
        gpk = fN1->GetNumberEntry()->GetIntNumber() -1;
        if (gpk>=(gMean.size()-1))    // keeping the "select" field limited to number of peaks found
        {
                fN1->SetNumber(gMean.size());
                gpk = fN1->GetNumberEntry()->GetIntNumber() -1;
        }
        if (!(gpk<gMean.size())) // just a precaution
        {
                gpk = gMean.size()-1;
        }
        TObject *des=gSearch->FindObject(selected);
        des->Delete();  // delete the last selection to create a new one
        float x[1] = {0};
        x[0] = gMean[gpk];
        float y[1] = {0};
        y[0] = gHeight[gpk];
        cout<<"Peak selected: "<< gMean[gpk]<<endl;
	gSearch->cd();
        selected= new TGraph(1,x,y);
        selected->SetMarkerStyle(23);
        selected->SetMarkerColor(1);
        selected->Draw("p");
        gSearch->Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Deletes the selected gamma peak and highlights it with a red triangular marker 
void gammaSearch::deletePeak()
{
	gSearch->cd();
        if (gMean.size()==1) {cout<<"Cannot delete the only Z peak."<<endl;}
        else
        {
		gpk = fN1->GetNumberEntry()->GetIntNumber() -1;
                float delx[1] = {gMean[gpk]};    // record deleted coordinates in temporary memory
                float dely[1] = {gHeight[gpk]};
                gMean.erase(gMean.begin()+gpk);  // delete the coordinates from the vectors that hold all peak info
                gHeight.erase(gHeight.begin()+gpk);
                cout<<"Erased"<<endl;
                TGraph *del = new TGraph(1,delx,dely); // plot the deleted peak as a red triangle
                del->SetMarkerStyle(23);
                del->SetMarkerColor(2);
                del->Draw("p");
                TObject *des=gSearch->FindObject(selected);   // delete the last selection to select the new peak
                des->Delete();
		// select the peak after the deleted one
                if (gpk>=(gMean.size()-1))
                {
                        fN1->SetNumber(gMean.size());
                        gpk = fN1->GetNumberEntry()->GetIntNumber() -1;
                }
                if (!(gpk<gMean.size()))
                {
                        gpk = gMean.size()-1;
                }
                float x[1] = {0};
                x[0] = gMean[gpk];
                float y[1] = {0};
                y[0] = gHeight[gpk];
                cout<<"Peak selected: "<< gMean[gpk]<<endl;
                selected= new TGraph(1,x,y);
                selected->SetMarkerStyle(23);
                selected->SetMarkerColor(1);
                selected->Draw("p");
                gSearch->Update();
        }
}

////////////////////////////////////////////////////////////////////////////////
/// This function is a wrapper function for performing the correlation of the found gamma peaks to those from the literature.
/// First, it calls the function that performs the correlations.
/// Then, it checks whether enough peaks were identified. It used the limit value from the limit vector - this is the percentage of the literature peaks that the program should identify for the correlation.
/// If not enough peaks were identified, it deletes the first found peak and does the correlation again, this time finding the ratios of the second found peak to all the next peaks.
/// The program will go through the peaks until it identifies the requested number of peaks and will label them with their detector channel and corresponding energy in keV.
/// If the program goes through all the peaks and doesn't identify enough peaks - an error message will appear. Most likely, a wrong calibration source was selected on the input page. It is also possible that the measurements are not accurate
/// If the program succeeds to find enough peaks, it will proceed to making a correlation plot of detector channels to energy
// CANNOT RESCALE MULTIPLE TIMES IN A ROW. NEED TO REDO THE SEARCH. OTHERWISE YOU'LL RESCALE THE ALREADY RESCALED HISTOGRAM
void gammaSearch::correlatePeaks()  // correlating channels to energies by looking at ratios between the first and other peaks
{
        corrected = (TH1F*)safeCopy->Clone();
	corrected->SetDirectory(0);
	found.clear();
	lit.clear();
	height.clear();
	ratios();
	while (found.size()<(allEnergy[iso].size()*limits[iso])) // more than limits% of peaks from the literature should be identified. If not enough were indentified, then the first found peak is not the first peak from the literature. Then the first peak will be deleted and the identification will start over 
	{
		found.clear(); // clearing all vectors so that no previous results are carried over to the next iteration of ratio search
		lit.clear();
		height.clear();
		gpk=0;
		deletePeak();     // perform the ratio search again
		if (gMean.size()==1)  // this is the last peak - need to break the cycle
		{
			cout<<"Something went wrong. You might have selected a wrong isotope, or you might need to increase the sensitivity during the peak search"<<endl;
			return;
		}
		ratios();
	}
	if (found.size()>=(allEnergy[iso].size()*limits[iso]))   // checking if the we found enough good ratios between peaks
	{
		for (unsigned int k=0; k<found.size();k++)  // labelling all peaks that passed the ratio test
		{
			gSearch->cd(1);
			text.push_back(new TLatex(found[k],height[k],Form("%ich=%ikeV",found[k],lit[k])));
	       		text[k]->SetTextSize(0.025);
 			text[k]->SetTextAngle(30.);
			text[k]->Draw();
        		gSearch->Update();
		}
		test = new TCanvas("E vs ch","E vs ch");   // making a plot of energy vs detector channels of the peaks that passed the ratio test
		TGraph *fit = new TGraph(lit.size(),found.data(),lit.data());
		fit->SetMarkerStyle(8);
		fit->SetTitle("Energy (keV) vs Channels");
		fit->GetXaxis()->SetTitle("Channel");
		fit->GetYaxis()->SetTitle("Energy (keV)");
		fit->Draw();
		TF1 *fitpk = new TF1("fitpk",calibration,20,4000,2); // 2 parameters  // fitting the calibration function
		fitpk->SetParLimits(0,0,1000000);
		fit->Fit("fitpk","p","Integral",found[0],found[found.size()-1]);
		//get parameters		
		m=fitpk->GetParameter(0);
		b=fitpk->GetParameter(1);
		TLatex function;  // displaying the fitted function equation on the canvas
		function.SetTextSize(0.025);
		function.SetTextAngle(0.);
		function.DrawLatex(found[0],lit[2],Form("E (keV) = %f * ch + %f",m,b));
		new ratioPeaks(gClient->GetRoot(),200,200);
	}
	else {cout<<"Something went wrong. Try manually deleting some false peaks and try again."<<endl;}
}

////////////////////////////////////////////////////////////////////////////////
/// This function performs the correlations between the found gamma peaks and those from the literature.
/// It looks at the ratios of the first gamma peak to all the other peaks in both the literature and the list of found peaks.
/// If the ratio between found peaks is within 1% of the ratio between the peaks from the literature, they are recorded in three vectors: 
///    - found peaks are recorded in the found vector
///    - peaks from the literature are recorded in the lit vector
///    - heights of the found peaks that passed the ratio test are recorded in the height vector 
void gammaSearch::ratios()
{
	lit.push_back(allEnergy[iso][0]);
	found.push_back(gMean[0]);
	height.push_back(gHeight[0]);
	for (unsigned int k=1; k<allEnergy[iso].size();k++)
	{
		float ratio1 = allEnergy[iso][0]/allEnergy[iso][k];   // ratio of the first literature energy peak to the rest 
		for (unsigned int m=1;m<gMean.size();m++)
		{
			float ratio2 = gMean[0]/gMean[m];   		// ratio of the first found peak to the rest 
			if ((ratio2-0.01*ratio2)< ratio1 && ratio1<(ratio2+0.01*ratio2))   // checking if the ratios are within 1% of each other. If yes - they pass the test
			{    
       				int a = round(gMean[m]);
       				int b = round(allEnergy[iso][k]);
				int c = round(gHeight[m]);
				found.push_back(a);
				lit.push_back(b);
				height.push_back(c);
			}
		//SUGGESTION FOR IMPROVEMENT: if two found peaks were correlated with the same literatre energy - delete both peaks
		}
	}
	cout<<"size"<<lit.size()<<","<<found.size()<<endl;
}

////////////////////////////////////////////////////////////////////////////////
/// This function makes a correlation plot of detector channels to energy and it marks just the peaks that have passed the correlation test on the histogram with green triangualar markers.
/// Here, the user can delete peaks that they don't want to fit. 
/// They can also go back to the previous window to delete the first peak if the correlation used wrong peaks. When the first peak is deleted and the correlation is launched again, it should successfully correlate the peaks.
ratioPeaks::ratioPeaks(const TGWindow *p, UInt_t w, UInt_t h)
{
        // create a main frame
        fMain1 = new TGMainFrame (p,w,h);
	
	// create a canvas widget
	fEcanvas1 = new TRootEmbeddedCanvas("Ecanvas1",fMain1, 800, 600);
	fMain1->AddFrame(fEcanvas1, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10,1));

        // create horizontal frames with text entry fields, labels and buttons
	TGHorizontalFrame *hframe2 = new TGHorizontalFrame(fMain1,200,40);
	fN0 = new TGNumberEntry(hframe2, 1,9,999,TGNumberFormat::kNESInteger,
                                                TGNumberFormat::kNEANonNegative,
                                                TGNumberFormat::kNELLimitMinMax,
                                                1,99999);
        fN0->Connect("ValueSet(Long_t)","ratioPeaks",this,"selectPeak()");
        fN0->GetNumberEntry()->Connect("ReturnPressed()","gammaSearch",this,"selectPeak()");
        TGLabel *fL0 = new TGLabel(hframe2,"Select the number of the peak to be removed (counted along the x-axis)");
        hframe2->AddFrame(fL0, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
        hframe2->AddFrame(fN0, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
	TGTextButton *del = new TGTextButton(hframe2,"&Delete This Peak");
        del->Connect("Clicked()","ratioPeaks",this,"deletePeak()");
        hframe2->AddFrame(del, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));

	fMain1->AddFrame(hframe2, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

        TGHorizontalFrame *hframe4 = new TGHorizontalFrame(fMain1,200,40);

        TGTextButton *gRedo = new TGTextButton(hframe4, "&Fit Peaks");
        gRedo->Connect("Clicked()","ratioPeaks",this,"fitPeaks()");
        hframe4->AddFrame(gRedo, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));

        TGTextButton *back = new TGTextButton(hframe4, "&Back");
        back->Connect("Clicked()","ratioPeaks",this,"back()");
        hframe4->AddFrame(back, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));

        TGTextButton *exit = new TGTextButton(hframe4, "&Exit", "gApplication->Terminate(0)");
        hframe4->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));
        fMain1->AddFrame(hframe4, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

        fMain1->SetCleanup(kDeepCleanup);

        // Set a name to the main frame
        fMain1->SetWindowName("Identified Peaks");
        // next 3 lines make the widgets visible     
        // Map all subwindows of main frame
        fMain1->MapSubwindows();
        // Initialize the layout algorithm
        fMain1->Resize(fMain1->GetDefaultSize());  // execute all layout specs for widgets before the top-level window is shown on screen
        // map main frame
        fMain1->MapWindow();
	
	for (unsigned int k=0;k<found.size();k++)
	{
		found[k]=m*found[k]+b;   // scaling the found peaks
	}
	auto ScaleX = [=](Double_t x){return (m*x+b);};
	keV=fEcanvas1->GetCanvas();
	keV->SetLogy();
	corrected->Draw();
	corrected->SetTitle("Gamma Energies");
	ScaleAxis(corrected->GetXaxis(), ScaleX);
	//plotting a scatter plot of all peaks
	gscat = new TGraph(found.size(),&found[0],&height[0]);
	gscat->SetMarkerStyle(23);
	gscat->SetMarkerColor(3);
	gscat->Draw("p");
//	for (int p=0;p<gfound;p++) {cout<<gMean[p]<<endl;}
	//creating a "selected" peak
	float x1[1]={0};
	x1[0]=found[0];
	float y1[1]={0};
	y1[0]=height[0];
	selected = new TGraph(1,x1,y1);
	selected->SetMarkerStyle(23);
	selected->SetMarkerColor(1);
	selected->Draw("p");
	keV->Update();
}

unsigned int epk =0; ///< This value stores hte index of a selected peak

////////////////////////////////////////////////////////////////////////////////
/// Selects a gamma peak and highlights it with a black triangle marker
void ratioPeaks::selectPeak()
{
        // get the number entry and convert it to 0 = 1st peak base
        TObject *sender = (TObject *)gTQSender;
        if (sender && sender->InheritsFrom("TGNumberEntry"))
        epk = fN0->GetNumberEntry()->GetIntNumber() -1;
        if (epk>=(found.size()-1))
        {
                fN0->SetNumber(found.size());
                epk = fN0->GetNumberEntry()->GetIntNumber() -1;
        }
        if (!(epk<found.size()))
        {
                epk = found.size()-1;
        }
        TObject *des=keV->FindObject(selected);
        des->Delete();  // delete the last selection to create a new one
        float x[1] = {0};
        x[0] = found[epk];
        float y[1] = {0};
        y[0] = height[epk];
        cout<<"Peak selected: "<< found[epk]<<endl;
	keV->Update();
	keV->cd();
        selected= new TGraph(1,x,y);
        selected->SetMarkerStyle(23);
        selected->SetMarkerColor(1);
        selected->Draw("p");
        keV->Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Deletes the selected gamma peak and highlights it with a red triangle marker
void ratioPeaks::deletePeak()
{
	keV->cd();
        if (found.size()==1) {cout<<"Cannot delete the only Z peak."<<endl;}
        else
        {
		epk = fN0->GetNumberEntry()->GetIntNumber() -1;
                int delx[1] = {found[epk]};    // record deleted coordinates in temporary memory
                int dely[1] = {height[epk]};
                found.erase(found.begin()+epk);  // delete the coordinates from the vectors that hold all peak info
                height.erase(height.begin()+epk);
		lit.erase(lit.begin()+epk);
                cout<<"Erased"<<endl;
                TGraph *del = new TGraph(1,delx,dely); // plot the deleted peak as a red triangle
                del->SetMarkerStyle(23);
                del->SetMarkerColor(2);
                del->Draw("p");
                TObject *des=keV->FindObject(selected);   // delete the last selection to select the new peak
                des->Delete();
                if (epk>=(found.size()-1))
                {
                        fN0->SetNumber(found.size());
                        epk = fN0->GetNumberEntry()->GetIntNumber() -1;
                }
                if (!(epk<found.size()))
                {
                        epk = found.size()-1;
                }
		// selecting the next peak
                float x[1] = {0};
                x[0] = found[epk];
                float y[1] = {0};
                y[0] = height[epk];
                cout<<"Peak selected: "<< found[epk]<<endl;
                selected= new TGraph(1,x,y);
                selected->SetMarkerStyle(23);
                selected->SetMarkerColor(1);
                selected->Draw("p");
                keV->Update();
        }
}

////////////////////////////////////////////////////////////////////////////////
/// This function is used to go back to before correlating the channels to energies. It closes the correlation plot and window, and goes back to the Gamma Search window.
/// This is used when the program identified wrong peaks. By going back, deleting the first peak and doing the correlation again, the user can fix this issue.
/// All the vectors that stored the values of the peaks that passed the calibration will be cleared when the correlation is redone.
void ratioPeaks::back()
{
	gSearch->cd();
	int index = 0;
	while (gSearch->FindObject(text[index])) 
	{cout<<"Going Back"<<endl;gSearch->Modified();gSearch->FindObject(text[index])->Delete(); gSearch->Update();index +=1;}
	if (test) {test->Close(); gSystem->ProcessEvents();}  // close the plot of energy vs the detector channel
	fMain1->CloseWindow();				      // close the window
	index = 0;
	text.clear();            // clear the labels of peaks that passed the ratio test
}

////////////////////////////////////////////////////////////////////////////////
/// This function is used to initiate the fitting of the gamma peaks that were not deleted. It clears all the used vectors to make sure that if the user went back, the old results are not carried over. 
void ratioPeaks::fitPeaks()
{
	energy.clear();
	Energy.clear();
	dEnergy.clear();
	Area.clear();
	dArea.clear();
	dEnergy.clear();
	yield.clear();
	dYield.clear();
	eff.clear();
	dEff.clear();
	new gammaFits(gClient->GetRoot(),200,200);   // initiate gamma fits
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a window to fit the gamma peaks and fits the first peak
/// Displays the energy at which the peak is located and the area of the peak (number of gamma particles detected in this peak)
gammaFits::gammaFits(const TGWindow *p, UInt_t w, UInt_t h)
{
        // create a main frame
        fMain = new TGMainFrame (p,w,h);
	
	// create a canvas widget
	fEcanvas2 = new TRootEmbeddedCanvas("Ecanvas2",fMain, 800, 600);
	fMain->AddFrame(fEcanvas2, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10,1));

        // create horizontal frames with text entry fields, labels and buttons
	TGHorizontalFrame *hframe1 = new TGHorizontalFrame(fMain,200,40);
	label = new TGLabel(hframe1,"Peak at Energy = ____________. Area = ______________");
	hframe1->AddFrame(label,new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
	fMain->AddFrame(hframe1, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

	TGHorizontalFrame *hframe2 = new TGHorizontalFrame(fMain,200,40);
        TGLabel *fL0 = new TGLabel(hframe2,"Use this peak for efficiency callibration? If no, the area can be entered in the entry field manually and click 'No'");
        hframe2->AddFrame(fL0, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
        fMain->AddFrame(hframe2, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

	TGHorizontalFrame *hframe3 = new TGHorizontalFrame(fMain,200,40);
	TGTextButton *yes = new TGTextButton(hframe3,"&Yes");
        yes->Connect("Clicked()","gammaFits",this,"gyes()");
        hframe3->AddFrame(yes, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));
	TGTextButton *no = new TGTextButton(hframe3,"&No");
        no->Connect("Clicked()","gammaFits",this,"gno()");
        hframe3->AddFrame(no, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));
	fN0= new TGNumberEntry(hframe3,0,9,100, TGNumberFormat::kNESRealTwo,
                                                 TGNumberFormat::kNEANonNegative,
                                                 TGNumberFormat::kNELLimitMinMax,
                                                 0,1000000);
	hframe3->AddFrame(fN0, new TGLayoutHints(kLHintsTop|kLHintsLeft,5,5,5,5));
	fMain->AddFrame(hframe3, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

        TGHorizontalFrame *hframe4 = new TGHorizontalFrame(fMain,200,40);
        TGTextButton *exit = new TGTextButton(hframe4, "&Exit", "gApplication->Terminate(0)");
        hframe4->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));
        fMain->AddFrame(hframe4, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

        fMain->SetCleanup(kDeepCleanup);

        // Set a name to the main frame
        fMain->SetWindowName("Gamma Fit");
        // next 3 lines make the widgets visible     
        // Map all subwindows of main frame
        fMain->MapSubwindows();
        // Initialize the layout algorithm
        fMain->Resize(fMain->GetDefaultSize());  // execute all layout specs for widgets before the top-level window is shown on screen
        // map main frame
        fMain->MapWindow();

	for (unsigned int o=0;o<found.size();o++) {cout<<found[o]<<endl; cout<<height[o]<<endl;}	
	// plot the histogram
	fitting = (TH1F*)safeCopy->Clone();
	auto ScaleX = [=](Double_t x){return (m*x+b);};
	gausFit=fEcanvas2->GetCanvas();
	gausFit->SetLogy();
	fitting->Draw();
	fitting->SetTitle("Gamma Energies");
	ScaleAxis(fitting->GetXaxis(), ScaleX);  
	// select the range of the fit
	float ll = found[peak_index]-14;
	float hl = found[peak_index]+14;
	fitting->GetXaxis()->SetRangeUser(ll,hl);
	// create the fit function
	TF1 *fitpeak = new TF1("fitpeak",gausbkg,ll,hl,6);  // 6 pars
	fitpeak->SetParLimits(0,0,1000000);
	fitpeak->SetParameter(0,height[peak_index]);   // height
	fitpeak->SetParameter(1,found[peak_index]);   // center for non skew
	fitpeak->SetParameter(2,1.5);   // st dev
//	fitpeak->SetParameter(3,0);   // skewness parameter
//	fitpeak->SetParameter(4,0);   // relative height of skewed gaussian
//	fitpeak->SetParameter(5,0.5);   // size of the step function
	fitpeak->SetParName(0,"Height");
	fitpeak->SetParName(1,"Center");
	fitpeak->SetParName(2,"Standard Deviation");
	fitpeak->SetParName(3,"Skewness Parameter");
	fitpeak->SetParName(4,"Relative Height of the Skewed Gaussian");
	fitpeak->SetParName(5,"Size of the Step Function");
	fitpeak->SetLineColor(2);
	fitpeak->Draw("same");
	
	// fit the function
	fitting->Fit("fitpeak",Form("%d",found[peak_index]),"Integral",ll,hl);
	//get parameters
	nrg = fitpeak->GetParameter(1);
	// create the background function to make it visible
	TF1 *fitbkg = new TF1("fitbkg",bkg,ll,hl,6);
	fitting->Fit("fitbkg",Form("%dbkg",found[peak_index]),"Integral",ll,hl);
	fitbkg->SetLineColor(1);
	fitbkg->Draw("same");
	area = 0;
	area = fitpeak->Integral(ll,hl)-fitbkg->Integral(ll,hl);  // calculate the integral of the peak to estimate the # of gamma particles in the peak
	gausFit->Update();
	double temp0 = (area-b)/m; // area scaled back to before channel calibration (needed bc otherwise the histogrma is stretched, and the area is stetched too)
	label->SetText(Form("Peak at Energy = %f. Area = %f",nrg,temp0));
}

////////////////////////////////////////////////////////////////////////////////
/// User approved the peak - add the area to the Area vector and move on to fit the next peak
void gammaFits::gyes()
{
	double temp0 = (area-b)/m; // area scaled back to before channel calibration
	Area.push_back(temp0);
	energy.push_back(nrg);
	cout<<"Written"<<nrg<<endl;
	gnext();       
}

////////////////////////////////////////////////////////////////////////////////
/// User disapproved the peak - check if the user entered a custom variable; if yes - use that variable, otherwise don't do anything and move onto the next peak
void gammaFits::gno()
{
	cout<<"Not written"<<nrg<<endl;
	float temp1 = ((fN0->GetNumberEntry()->GetNumber())-b)/m;
	if (!(round(temp1)==0))  // if the user entered a custom value instead
	{
		Area.push_back(temp1);
		energy.push_back(nrg); 
		cout<<"Manually entered area written "<<temp1<<endl;
		fN0->SetNumber(0);
	}
	gnext();
}

////////////////////////////////////////////////////////////////////////////////
/// Fit the next peak and update the label with a new peak energy and # of gamma particles. If all peaks were fitted, move onto fitting the efficiency function
void gammaFits::gnext()
{
	if (peak_index+1<found.size())
	{
		peak_index +=1;
		float ll = found[peak_index]-14;
		float hl = found[peak_index]+14;
		TF1 *fitpeak = new TF1("fitpeak",gausbkg,ll,hl,6);  // 6 pars
		fitpeak->SetParLimits(0,0,1000000);
		fitpeak->SetParameter(0,height[peak_index]);   // height
		fitpeak->SetParameter(1,found[peak_index]);   // center for non skew
		fitpeak->SetParameter(2,1.5);   // st dev
		fitpeak->SetParameter(3,10);   // skewness parameter
		fitpeak->SetParName(0,"Height");
		fitpeak->SetParName(1,"Center");
		fitpeak->SetParName(2,"Standard Deviation");
		fitpeak->SetParName(3,"Background, vertical shift");
		fitpeak->SetParName(4,"Background, horizontal stretch");
		fitpeak->SetParName(5,"Background, vertical stretch");
		fitpeak->SetLineColor(2);
		fitpeak->Draw("same");
		
		fitting->GetXaxis()->SetRangeUser(ll,hl);
		fitting->Fit("fitpeak",Form("%d",found[peak_index]),"Integral",ll,hl);
		//get parameters
		float p1 = fitpeak->GetParameter(1);
		float p3 = fitpeak->GetParameter(3);
		float p4 = fitpeak->GetParameter(4);
		float p5 = fitpeak->GetParameter(5);
		nrg = fitpeak->GetParameter(1);
		TF1 *fitbkg = new TF1("fitbkg",bkg,ll,hl,6);
		fitbkg->FixParameter(1,p1);
		fitbkg->FixParameter(3,p3);
		fitbkg->FixParameter(4,p4);
		fitbkg->FixParameter(5,p5);
		fitbkg->SetLineColor(1);
		fitbkg->Draw("same");
		area = 0;
		area = fitpeak->Integral(ll,hl)-fitbkg->Integral(ll,hl);
		gausFit->Update();
		double temp0 = (area-b)/m; // area scaled back to before channel calibration
		label->SetText(Form("Peak at Energy = %f. Area = %f",nrg,temp0));
	}
	else 
	{
		fitting->GetXaxis()->UnZoom();
		fitting->GetYaxis()->UnZoom();
		gausFit->Update();
		gPad->Modified();
		gPad->Update();
		label->SetText("All peaks were already plotted."); 
		label->SetTextColor(kRed);
		cout<<peak_index<<"All peaks were already plotted"<<endl;
		new efficiency(gClient->GetRoot(),200,200);
	}

}

////////////////////////////////////////////////////////////////////////////////
/// Creates a window to fit the efficiency function. Calls the function that will plot and fit the efficiency curve.
efficiency::efficiency(const TGWindow *p, UInt_t w, UInt_t h)
{
        // create a main frame
        fMain = new TGMainFrame (p,w,h);
	
	// create a canvas widget
	fEcanvas3 = new TRootEmbeddedCanvas("Ecanvas3",fMain, 800, 600);
	fMain->AddFrame(fEcanvas3, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10,1));

        // create horizontal frames with text entry fields, labels and buttons
	TGHorizontalFrame *hframe5 = new TGHorizontalFrame(fMain,200,40);
        TGTextButton *exit = new TGTextButton(hframe5, "&Exit", "gApplication->Terminate(0)");
        hframe5->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));
	fMain->AddFrame(hframe5, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
	
        fMain->SetCleanup(kDeepCleanup);

        // Set a name to the main frame
        fMain->SetWindowName("Efficiency");
        // next 3 lines make the widgets visible     
        // Map all subwindows of main frame
        fMain->MapSubwindows();
        // Initialize the layout algorithm
        fMain->Resize(fMain->GetDefaultSize());  // execute all layout specs for widgets before the top-level window is shown on screen
        // map main frame
        fMain->MapWindow();
	
	effi = fEcanvas3->GetCanvas();
	plot();
}

////////////////////////////////////////////////////////////////////////////////
/// Plots and fits the efficiency curve.
/// First, this function makes vectors for energy, yield and their uncertainties for the user-approved peaks. If the found peak is not within the 2.3 keV of the literature value, it will be ignored.
/// Then the efficiencies with uncertainties are calculated for each peak.
/// All efficiency points are plotted on the canvas. Then the effFunc is used to fit the efficiency curve. 
void efficiency::plot()
{
	for (unsigned int i=0;i<energy.size();i++){for (unsigned int j=0;j<allEnergy[iso].size();j++)
	{
		if ((energy[i]-2.3)<allEnergy[iso][j] && allEnergy[iso][j]<(energy[i]+2.3))   // making energy and yield vectors for all peaks that passed the ratio test
		{
			Energy.push_back(allEnergy[iso][j]);
			dEnergy.push_back(alldEnergy[iso][j]);
			yield.push_back(allYield[iso][j]);
			dYield.push_back(alldYield[iso][j]);
			break;
		}
	}}
	//calculating the efficiency
	for (unsigned int k=0;k<Energy.size();k++)
	{
cout<<"Area: "<<Area[k]<<" Act: " <<activity<< " Time: "<<Time<< "Yield: "<<yield[k]<<endl;
		eff.push_back(Area[k]/activity/Time/yield[k]);
		dArea.push_back(sqrt(Area[k]));
		dEff.push_back(eff[k]*sqrt(1/Area[k]+pow((dYield[k]/yield[k]),2)+pow((dtime/Time),2)+pow((dact/activity),2)));
	} 
	effi->cd();
	effi->SetLogy();
	effi->SetLogx();
	effi->Draw();
	// plotting the efficiency
	TGraphErrors *graph = new TGraphErrors(Energy.size(),&Energy[0],&eff[0],&dEnergy[0],&dEff[0]);
	graph->SetTitle("Efficiency vs Energy");
	graph->GetYaxis()->SetTitle("Efficiency");
	graph->GetXaxis()->SetTitle("Energy");
	graph->Draw();
	effi->Update();
	// making a fit function to fit the efficiency
	TF1 *fitplot = new TF1("fitplot", effFunc, 10,7000,6); // 6 pars   - CHANGE 6 IF NEW PARAMETERS WERE ADDED TO THE EFFICIENCY FUNCTION
	fitplot->SetParLimits(0,0,100000);
	fitplot->SetParameter(0,0);
	fitplot->SetParameter(1,0);
	fitplot->SetParameter(2,0);
	// fitting the efficiency 
	graph->Fit("fitplot","p","Integral",Energy[0],Energy[Energy.size()-1]);
	float p0 = fitplot->GetParameter(0);
	float p1 = fitplot->GetParameter(1);
	float p2 = fitplot->GetParameter(2);
	float p3 = fitplot->GetParameter(3);
	float p4 = fitplot->GetParameter(4);
	float p5 = fitplot->GetParameter(5);    /// USE THIS TEMPLATE FOR EVERY NEW PARAMETER : float p# = fitplot->GetParameter(#), # is the parameter number
	gStyle->SetOptFit();
	// displaying the formula of the fitted efficiency function on the screen
	TLatex formula;
	formula.SetTextSize(0.025);
	formula.SetTextAngle(0.);
	formula.DrawLatex((graph->GetXaxis()->GetXmin())*1.05,(*min_element(eff.begin(), eff.end()))*0.95,Form("Eff=exp((%f+%fln(E)+%f(ln(E))^{2})#frac{2}{#pi}atan(exp(%f+%fln(E)+%f(ln(E))^{2}))-25)",p0,p1,p2,p3,p4,p5));  // CHANGE THIS EXPRESSION IF THE EFFICIENCY FUNCTION WAS CHANGED
	effi->Update();
	for (unsigned int n=0;n<eff.size();n++)
	{
		cout<<"Energy: "<<Energy[n]<<"Efficiency :"<<eff[n]<<endl;
	}
}

////////////////////////////////////////////////////////////////////////////////
/// Deconstructor for the input class.
input::~input()
{
        fMain->Cleanup();
        delete fMain;
}
////////////////////////////////////////////////////////////////////////////////
/// Deconstructor for the gammaSearch class.
gammaSearch::~gammaSearch()
{
	fMain->Cleanup();
	delete fMain;
}
////////////////////////////////////////////////////////////////////////////////
/// Deconstructor for the ratioPeaks class.
ratioPeaks::~ratioPeaks()
{
	fMain1->Cleanup();
	delete fMain1;
}
////////////////////////////////////////////////////////////////////////////////
/// Deconstructor for the gammaFits class.
gammaFits::~gammaFits()
{
	fMain->Cleanup();
	delete fMain;
}
////////////////////////////////////////////////////////////////////////////////
/// Deconstructor for the efficiency class.
efficiency::~efficiency()
{
	fMain->Cleanup();
	delete fMain;
}

////////////////////////////////////////////////////////////////////////////////
/// This is a linear function used to correlate the detector channels to energy levels.
Double_t calibration(Double_t *x, Double_t *par)
{
	Double_t func = par[0]*x[0]+par[1];
	return func;
}

////////////////////////////////////////////////////////////////////////////////
/// Gauss function used to fit gamma peaks together with bkg.
Double_t gausbkg(Double_t *x, Double_t *par)   // fitting function with curvy bkg
{
	// x[0]
	// par[0] = height
	// par[1] = center for non skew
	// par[2] = st dev of the gaus
	// par[3] = vertical shift of the background 
	// par[4] = horizontal stretch of the background function
	// par[5] = vertical stretch  of the background function
	
	Double_t gauss = par[0] *  TMath::Gaus(x[0], par[1], par[2]) +bkg(x,par);
	return gauss;
}

////////////////////////////////////////////////////////////////////////////////
/// Used by the gausbkg to fit background when fitting gamma peaks. Uses an error function for better estimate of the background "step" behaviour.
Double_t bkg(Double_t *x, Double_t *par)       // curvy bkg for the gaus function
{
	Double_t step_func =par[3]+ par[5] *TMath::Erfc((x[0]-par[1])/par[4]);
	return step_func;
}

////////////////////////////////////////////////////////////////////////////////
/// This is the efficiency function used to fit the efficiency curve.
Double_t effFunc(Double_t *x, Double_t *par)
{
	Double_t pi = TMath::Pi();
	Double_t z=TMath::Log(x[0]);
//	Double_t f = par[0]+par[1]*z+par[2]*z*z;  // simplified
	Double_t f = TMath::Exp((par[0]+par[1]*z+par[2]*z*z)*2/pi*(TMath::ATan(TMath::Exp(par[3]+par[4]*z+par[5]*z*z)))-25.);
	return f;
}

////////////////////////////////////////////////////////////////////////////////
/// Used to scale the X axis.
/// Reference: https://root-forum.cern.ch/
void ScaleAxis(TAxis *a, std::function<Double_t(Double_t)> Scale)
{
	if (!a) return; // just a precaution
	if (a->GetXbins()->GetSize())
	{
		TArrayD X(*(a->GetXbins()));
		for(Int_t i = 0; i < X.GetSize(); i++) X[i] = Scale(X[i]);
		a->Set((X.GetSize() - 1), X.GetArray()); // new Xbins
	}
	else
	{
		a->Set(a->GetNbins(),
			Scale(a->GetXmin()), // new Xmin
			Scale(a->GetXmax()) ); // new Xmax
	}
	return;
}

////////////////////////////////////////////////////////////////////////////////
/// Initiates the program. IMPORTANT: DO NOT CHANGE THIS FUNCTION.
int main(int argc, char *argv[])
{
	TApplication theApp("App", &argc, argv);
	std::cout << "Launching... " << std::endl;
	new input(gClient->GetRoot(),200,200);
	theApp.Run();
	return 0;
}