////////////////////////////////////////////////
//Author: Sophia Devinyak (s.devinyak@gmail.com , sdevinya@uwaterloo.ca)
//Date created: November 2021 at TRIUMF under the supervision of Dr. Roger Caballero-Folch

READ ME


BEFORE USAGE
This program can only be run through the terminal. It cannot be launched by double clicking. root6 must be started before launching the program.

Whenever EfficiencyCalibrator is downloaded for the first time or moved to a different directory, the user must first run the CompileGammaCalibration.sh file to make sure that there are no errors when using the program. Whenever the user sees the following error, please run the CompileGammaCalibration.sh script before using the EfficiencyCalibrator:

Error in cling::AutoLoadingVisitor::InsertIntoAutoLoadingState:
   Missing FileEntry for input.h
   requested to autoload type input

This script can be run by typing out ./CompileGammaCalibration.sh in the terminal.


PREPARING THE FILE
This program analyzes a 2D histogram, where x-axis is the channel of the detector, and y-axis is the number of detected counts. Save your TCanvas containing the TH1F* histogram as a .root file before using EfficiencyCalibrator. 

INPUTTING THE INFORMATION
By default, the uncertainty in time is 1 second. The default value will be used if 0 is entered in the number entry field. 

GAMMA PEAKS SEARCH
The sensitivity of the search will determine how many peaks is identivied. The sensitivity is approximately the height of the smallest peak divided by the height of the tallest peak (therefore, smaller number will detect smaller peaks). Any false peaks can be deleted manually if desired, but most of the times it is not required for calibrating detector channels to energy levels. By clicking "Correlate Found Peaks", the program will ignore the peaks that are noise. This function looks at the ratios between the found peaks and compares them to the ratios of the peaks from the literature to correlate them. If wrong peaks were detected, the user can go back, delete the first peak and try again. This should usually fix the correlation.

SELECTING PEAKS
All found peaks are marked with green triangles. The selected peak is marked with a black triangle. The user can move between peaks by using the arrows on the entry field or by typing in the number of the desired peak and clicking "Enter". The peak can be deleted from calibration. The deleted peak will be marked with a black triangle.
	NOTE: a certain amount of peaks must match with the peaks from the literature. If too many peaks were deleted, the program will not be able 		to make the channel->energy calibration.

GAMMA PEAKS FIT
When fitting the gamma peaks, the user can manually select which ones to use in the efficiency calculation by using the "Yes"/"No" buttons. After reaching the last peak, the program will move on to plotting and fitting the efficiency. If the user is not satisfied, the area can be entered manually in the entry field. When the user clicks "No", if the entry field is 0, then the peak will be ignored. If the entry field contains user-entered value, then this value will be used instead.

The area to enter manually should be extracted from the calibrated spectrum. To do so, the user MUST follow these instructions, or the calibration will not work.
	i. Save the calibrated spectrum by right clicking outside the frame on the canvas and selecting SaveAs. 
	ii. Enter the file name.root in the first field. The second field can be left blank.
	iii. Open the saved file in a separate terminal. Fit a Gaussian to the peak, find its name by right clicking on the curve.
	iv. The area of the peak can be calculated by using the command *name of the curve*->Integral(lower limit of the x-axis, higher limit), then fitting a first or second degree polynomial to the background, finding its integral in the same limits, and subtracting this value from the Gauss curve integral.

CORRELATION FUNCTION FITTING
Make sure that you do not delete important datapoints for fitting to ensure that the fit works properly. Make sure that you have a couple of data points in the lower energy levels, a couple in higher levels, and some in the middle.

ADDING NEW CALIBRATION SOURCES
To add new calibration sources, you need to edit the DetectorEfficiency.C file. At the beginning of the script, you will find the vectors containing the information about the calibration sources. Append the vectors allIsotopes, allHalfLives, dallHalfLives, limits, allEnergy, alldEnergy allYield, alldYield with information about the isotope. When appending allEnergy, alldEnergy allYield and alldYield, erase the last };, add a comma, then in a new line type in { _____ }}; replacing _____ with respective values.

FITTING FUNCTION
If you wish to modify the function used to fit the efficiency curve, it can be done by editing the effFunc function at the end of the DetectorEfficiency.C file. Then go to the definition of the function efficiency::plot() and follow the instructions in the comments that were written in all caps. After editing, make sure to run the CompileGammaCalibration.sh script before using the program to implement the changes.

