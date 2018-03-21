#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>

#define THRESHOLDSTART 3
#define THRESHOLDEND 11.5
#define THRESHOLDSET false
#define THRESHOLDT2 362

//perform simulated real-time ripple detection on this file
#define DATAINPUT "/home/shayok/Documents/Code/RippleDetectionAnalysis/Cavaradossi/paperData/syntheticRippleAnalysis/smoothed_envelope_simulated.out"

//Since we're doing a simulated signal we're only gonna have 1 bootstraps
//all of the signal looks similar so it'll be representative if we just look at 
//the data all together.
#define MINUTECOUNT 15 //how many minutes are the bootstrap samples
#define BOOTSTRAPFILENAME "/home/shayok/Documents/Code/RippleDetectionAnalysis/Cavaradossi/paperData/syntheticRippleAnalysis/bootstrapStartPoints.out"

#define RIPPLESTARTBOUND "/home/shayok/Documents/Code/RippleDetectionAnalysis/Cavaradossi/paperData/syntheticRippleAnalysis/rippleBoundsStart.out"
#define RIPPLEENDBOUND "/home/shayok/Documents/Code/RippleDetectionAnalysis/Cavaradossi/paperData/syntheticRippleAnalysis/rippleBoundsEnd.out"


double calcMean(std::vector<double> arrrayForEst)
{
    double sums = 0;
    for(unsigned int x = 0; x<arrrayForEst.size();++x){
        sums += arrrayForEst[x];
    }
    double mean = sums/arrrayForEst.size();
    return mean;
}

double calcSTD(std::vector<double> arrrayForEst)
{
    double mean = calcMean(arrrayForEst);
    double sums = 0;
    for(unsigned int x = 0; x<arrrayForEst.size(); ++x)
        sums += (arrrayForEst[x]-mean) * (arrrayForEst[x]-mean);

    sums = sums/arrrayForEst.size();
    double standardDeviation = std::sqrt(sums);
    return standardDeviation;
}

/**
 * @brief calcTPRate
 * @param rippleStart
 * @param rippleEnd
 * @param detections
 * @param startIndex
 * @param endIndex
 * @param blockLengthIndex
 * @return
 * loop through all putative ripples (within startIndex and endIndex) and then loop through all detected ripples making sure and count how many are within the ripple bounds (ignore
 * ripple bounds that are within 200 ms of the previous one and count detections that are within 200ms of the start of the ripple bound)
 */
double calcTPRate(std::vector<int> rippleStart, std::vector<int> rippleEnd, std::vector<int> detections, int startIndex, int endIndex, int blockLengthIndex)
{
    //determine start and end index of putative ripple detections
    int i=0;
    if(i == rippleStart.size()){
        std::cout<<"Error no ripples within the 15-20 minutes?"<<std::endl;
        return -1;
    }
    while(rippleEnd[i]<startIndex && i<rippleStart.size())//use rippleEnd here to account for cases where the startIndex is in the middle of a ripple!
        ++i;
    int end = i;
    while(rippleEnd[end]<endIndex && end < rippleEnd.size())
        ++end;

    //determine start index of actual ripple detections
    int j=0;
    while(detections[j]<rippleStart[i] && j<detections.size())
        ++j;
    if(j == detections.size()){
        std::cout<<"Error no ripples detections within the 15-20 minutes?"<<std::endl;
        return -1;
    }

    int totalRipples = end - i;
    int tpCount = 0;
    bool firsttime = false;
    for(i; i<end; ++i){
        while(detections[j]<rippleStart[i]){ //I don't care about multiple detections per ripple block since they'd be suppressed if we were actually stimulating
            ++j;
            if(j > detections.size()){
                j = detections.size()-1;
                break;
            }
        }
        //let's start by discounting any putative ripple bounds that are within 200 ms of the previous one
        if(rippleStart[i] - rippleEnd[i-1] < blockLengthIndex && !(detections[j]>rippleStart[i] && detections[j]<rippleEnd[i])){
            --totalRipples;
        }
        //on time detections...best case
        else if((rippleStart[i] < detections[j]) && (rippleEnd[i] > detections[j])){
            ++tpCount;
        }
        if(j+1 == detections.size())
            break;
    }
    return (double)tpCount/totalRipples ;
}

/**
 * @brief calcFPRate
 * @param rippleStart
 * @param rippleEnd
 * @param detections
 * @param startIndex
 * @param endIndex
 * @param blockLengthIndex
 * @param time
 * @return false stimulation rate per unit of time parameter passed (e.g. if time is 12 minutes then it'll be per minute. if time is 12 seconds then it'll be per second.)
 */
double calcFPRate(std::vector<int> rippleStart, std::vector<int> rippleEnd, std::vector<int> detections, int startIndex, int endIndex, int blockLengthIndex, double time)
{
    //determine start index of putative ripple detections
    int i=0;
    while(rippleEnd[i]<startIndex && i<rippleStart.size())//use rippleEnd here to account for cases where the startIndex is in the middle of a ripple!
        ++i;
    if(i == rippleStart.size()){
        std::cout<<"Error no ripples within the 20 minutes?"<<std::endl;
        return -1;
    }

    int enddd = i;
    while(rippleEnd[enddd]<endIndex && enddd<rippleEnd.size())
        ++enddd;

    //determine start and end index of actual ripple detections
    int j=0;
    while(detections[j]<rippleEnd[i] && j<detections.size())
        ++j;
    if(j == detections.size()){
        std::cout<<"Error no ripples detections within the 20 minutes?"<<std::endl;
        return -1;
    }
    while(detections[j]>rippleEnd[i])
        --j;
    ++j;

    int STAAAARRTR = j;
    int end = j;
    while(detections[end]<endIndex && end < detections.size())
        ++end;


    int totalDetections = end - j;
    int fpCount = 0;
    for(i;i<enddd;++i){
        while (rippleStart[i] < detections[j] && detections[j] < rippleEnd[i]) {
            ++j;
        }
        if(i!=enddd){
            while(detections[j]<rippleStart[i+1] && j<end){
                ++j;
                ++fpCount;
                if(rippleStart[i+1] - detections[j] < blockLengthIndex && rippleStart[i+1]>detections[j]){ //account for early detections
                    --fpCount;
                }
            }
        }
    }
    return (double)fpCount/(time);
}

/**
 * @brief calcFPRate
 * @param rippleStart
 * @param rippleEnd
 * @param detections
 * @param startIndex
 * @param endIndex
 * @param blockLengthIndex
 * @return false positive percentage out of total online detections within time indecies
 */
double calcFPPercentage(std::vector<int> rippleStart, std::vector<int> rippleEnd, std::vector<int> detections, int startIndex, int endIndex, int blockLengthIndex)
{
    //determine start index of putative ripple detections
    int i=0;
    while(rippleEnd[i]<startIndex && i<rippleStart.size())//use rippleEnd here to account for cases where the startIndex is in the middle of a ripple!
        ++i;
    if(i == rippleStart.size()){
        std::cout<<"Error no ripples within the 20 minutes?"<<std::endl;
        return -1;
    }

    int enddd = i;
    while(rippleEnd[enddd]<endIndex && enddd<rippleEnd.size())
        ++enddd;

    //determine start and end index of actual ripple detections
    int j=0;
    while(detections[j]<rippleEnd[i] && j<detections.size())
        ++j;
    if(j == detections.size()){
        std::cout<<"Error no ripples detections within the 20 minutes?"<<std::endl;
        return -1;
    }
    while(detections[j]>rippleEnd[i])
        --j;
    ++j;

    int STAAAARRTR = j;
    int end = j;
    while(detections[end]<endIndex && end < detections.size())
        ++end;


    int totalDetections = end - j - 1;
    int fpCount = 0;
    for(i;i<enddd;++i){
        while (rippleStart[i] < detections[j] && detections[j] < rippleEnd[i]) {
            ++j;
        }
        if(i!=enddd){
            while(detections[j]<rippleStart[i+1] && j<end){
                ++j;
                ++fpCount;
                if(rippleStart[i+1] - detections[j] < blockLengthIndex && rippleStart[i+1]>detections[j]){ //account for early detections
                    --fpCount;
                }
            }
        }
    }
    return (double)fpCount/totalDetections ;
}

/**
 * @brief calcDetectLatencies
 * @param rippleStart
 * @param rippleEnd
 * @param detections
 * @param startIndex
 * @param endIndex
 * @param blockLengthIndex
 * @return average detection latency in ms
 */
double calcDetectLatencies(std::vector<int> rippleStart, std::vector<int> rippleEnd, std::vector<int> detections, int startIndex, int endIndex, int blockLengthIndex)
{
    //determine start and end index of putative ripple detections
    int i=0;
    while(rippleEnd[i]<startIndex && i<rippleStart.size())//use rippleEnd here to account for cases where the startIndex is in the middle of a ripple!
        ++i;
    if(i == rippleStart.size()){
        std::cout<<"Error no ripples within the 20 minutes?"<<std::endl;
        return -1;
    }
    int end = i;
    while(rippleEnd[end]<endIndex && end < rippleEnd.size())
        ++end;

    //determine start index of actual ripple detections
    int j=0;
    while(detections[j]<rippleStart[i] && j<detections.size())
        ++j;
    if(j == detections.size()){
        std::cout<<"Error no ripples detections within the 20 minutes?"<<std::endl;
        return -1;
    }

    std::vector<double> detectionLatencies;//detection latencies in ms

    for(i; i<end; ++i){
        if(i>0){
            if(rippleStart[i]-rippleEnd[i-1] > blockLengthIndex && (detections[j]>rippleStart[i] && detections[j]<rippleEnd[i])){
                if((detections[j]>rippleStart[i] && detections[j]<rippleEnd[i])){
                    detectionLatencies.push_back((detections[j]-rippleStart[i])/3);
                }
            }
        }
        else{
            if((detections[j]>rippleStart[i] && detections[j]<rippleEnd[i])){
                detectionLatencies.push_back((detections[j]-rippleStart[i])/3);
            }
        }
        while(detections[j]<rippleStart[i+1]){
            ++j;
        }
    }

    return calcMean(detectionLatencies) ;
}

/**
 * @brief calcDetectLatencies
 * @param rippleStart
 * @param rippleEnd
 * @param detections
 * @param startIndex
 * @param endIndex
 * @param blockLengthIndex
 * @return average relative detection latency: percent of event that has transpired before ripple was detected
 */
double calcRelativeDetectLatencies(std::vector<int> rippleStart, std::vector<int> rippleEnd, std::vector<int> detections, int startIndex, int endIndex, int blockLengthIndex)
{
    //determine start and end index of putative ripple detections
    int i=0;
    while(rippleEnd[i]<startIndex && i<rippleStart.size())//use rippleEnd here to account for cases where the startIndex is in the middle of a ripple!
        ++i;
    if(i == rippleStart.size()){
        std::cout<<"Error no ripples within the 20 minutes?"<<std::endl;
        return -1;
    }
    int end = i;
    while(rippleEnd[end]<endIndex && end < rippleEnd.size())
        ++end;

    //determine start index of actual ripple detections
    int j=0;
    while(detections[j]<rippleStart[i] && j<detections.size())
        ++j;
    if(j == detections.size()){
        std::cout<<"Error no ripples detections within the 20 minutes?"<<std::endl;
        return -1;
    }

    std::vector<double> detectionLatencies;//detection latencies in ms

    for(i; i<end; ++i){
        if(i>0){
            if(rippleStart[i]-rippleEnd[i-1] > blockLengthIndex && (detections[j]>rippleStart[i] && detections[j]<rippleEnd[i])){
                if((detections[j]>rippleStart[i] && detections[j]<rippleEnd[i])){
                    double moo = detections[j]-rippleStart[i];
                    double cow = rippleEnd[i]-rippleStart[i];
                    detectionLatencies.push_back(moo*100/cow);
                }
            }
        }
        else{
            if((detections[j]>rippleStart[i] && detections[j]<rippleEnd[i])){
                double cow = rippleEnd[i]-rippleStart[i];
                double moo = detections[j]-rippleStart[i];
                detectionLatencies.push_back(moo*100/cow);
            }
        }
        while(detections[j]<rippleStart[i+1]){
            ++j;
        }
    }
    return calcMean(detectionLatencies) ;
}

int main(int argc, char *argv[])
{
    std::string SIMDETECTIONFILENAME = argv[1];
    SIMDETECTIONFILENAME += "/simDetectionsSingleChan";

    std::string TPRATEFILENAME = argv[1];
    TPRATEFILENAME += "/tpRate";
    
    std::string FPRATEFILENAME = argv[1];
    FPRATEFILENAME += "/fpRate";
    
    std::string FPPERCENTAGEFILENAME = argv[1];
    FPPERCENTAGEFILENAME += "/fpPercent";
    
    std::string DETECTIONLATENCYFILENAME = argv[1];
    DETECTIONLATENCYFILENAME += "/detectionLatency";
    
    std::string RELATIVEDETECTIONLATENCYFILENAME = argv[1];
    RELATIVEDETECTIONLATENCYFILENAME += "/relativeDetectionLatency";

    // int minIndexCount = 3000*60*MINUTECOUNT; //mincount*60*FS...FS=3kHz
    int minIndexCount = 14408422; //length of dataset
    /*
    load up bootstrap start times...this is kinda useless but I just copied 
    and pasted to keep it consistent with the in vivo analysis. useless
    because bootstrap start time is 0. one sample. 
    */
    std::ifstream startIntervals(BOOTSTRAPFILENAME); 
    std::vector <int> bootstrapStartIntervals;
    std::string line;
    if(startIntervals.is_open()){
        while(std::getline(startIntervals,line)){
            bootstrapStartIntervals.push_back(std::stod(line));
        }
    }
    else{
        std::cout << "Error opening file time intervals"<<std::endl;
        return 1;
    }

    //load up ripple onset times
    std::ifstream startBounds(RIPPLESTARTBOUND);
    std::vector <int> rippleBoundStart;
    if(startBounds.is_open()){
        while(std::getline(startBounds,line)){
            rippleBoundStart.push_back(std::stod(line));
        }
    }
    else{
        std::cout << "Error opening file ripple bounds start"<<std::endl;
        return 1;
    }

    //load up ripple end bounds
    std::ifstream endBounds(RIPPLEENDBOUND);
    std::vector <int> rippleBoundEnd;
    if(endBounds.is_open()){
        while(std::getline(endBounds,line)){
            rippleBoundEnd.push_back(std::stod(line));
        }
    }
    else{
        std::cout << "Error opening file ripple bounds end"<<std::endl;
        return 1;
    }

    int blockLength = 600; //200ms block after detection or 600 indexes at sampling rate 3kHz

    for(double x = THRESHOLDSTART; x<THRESHOLDEND; x+=0.25){ //loop through a bunch of thresholds, perform ripple detections, and evaluate
        std::cout << x<<std::endl;
        std::ostringstream strs;
        strs << x*100;

        std::ofstream myfile;
        std::string fileNameee = SIMDETECTIONFILENAME;

        
        fileNameee += strs.str() + "out";
        myfile.open(fileNameee, std::ofstream::out | std::ofstream::app);

        //Read in Files

        std::ifstream dataInT2(DATAINPUT);

        std::vector <double> smoothed_envelopeT2;
        std::vector <int> detectionTimeIndexes;

        double thresholdT2;

        if(dataInT2.is_open()){
            while(std::getline(dataInT2,line)){
                smoothed_envelopeT2.push_back(std::stod(line));
            }
            if(THRESHOLDSET){
                thresholdT2 = THRESHOLDT2;
            }
            else{
                // thresholdT2 = calcMean(smoothed_envelopeT2)+(x*calcSTD(smoothed_envelopeT2));
                thresholdT2 = (7.53744650311023)+(x*3.5462055039255231);
                std::cout << thresholdT2<<std::endl;
            }
        }
        else{
            std::cout <<"Error opening file T2"<<std::endl;
            return 1;
        }

        //Hunt for ripples on two channels
        int i=0;
        while (i<smoothed_envelopeT2.size()){ //loop through all elements of channel
            if(smoothed_envelopeT2[i] > thresholdT2){
                myfile << i << " " << i+10 <<"\n";
                myfile.flush();
                detectionTimeIndexes.push_back(i);
                i+=blockLength;
            }
            else{
                ++i;
            }
        }

        std::ofstream tpRate;
        fileNameee = TPRATEFILENAME;

        fileNameee += strs.str() + "out";
        tpRate.open(fileNameee, std::ofstream::out | std::ofstream::app);

        std::ofstream fpRate;
        fileNameee = FPRATEFILENAME;

        fileNameee += strs.str() + "out";
        fpRate.open(fileNameee, std::ofstream::out | std::ofstream::app);

        std::ofstream detectionLatency;
        fileNameee =DETECTIONLATENCYFILENAME;

        fileNameee += strs.str() + "out";
        detectionLatency.open(fileNameee, std::ofstream::out | std::ofstream::app);

        std::ofstream relativeDetectionLatency;
        fileNameee =RELATIVEDETECTIONLATENCYFILENAME;

        fileNameee += strs.str() + "out";
        relativeDetectionLatency.open(fileNameee, std::ofstream::out | std::ofstream::app);

        std::ofstream fpPercentage;
        fileNameee = FPPERCENTAGEFILENAME;

        fileNameee += strs.str() + "out";
        fpPercentage.open(fileNameee, std::ofstream::out | std::ofstream::app);

        //loop through all bootstrap samples, evaluate ripple detrections and print results to file
        for(int xa = 0; xa<bootstrapStartIntervals.size();++xa){
            tpRate << calcTPRate(rippleBoundStart, rippleBoundEnd, detectionTimeIndexes, bootstrapStartIntervals[xa], bootstrapStartIntervals[xa]+minIndexCount, blockLength) <<'\n';
            fpRate << calcFPRate(rippleBoundStart, rippleBoundEnd, detectionTimeIndexes, bootstrapStartIntervals[xa], bootstrapStartIntervals[xa]+minIndexCount, blockLength,15) <<'\n';
            fpPercentage << calcFPPercentage(rippleBoundStart, rippleBoundEnd, detectionTimeIndexes, bootstrapStartIntervals[xa], bootstrapStartIntervals[xa]+minIndexCount, blockLength) <<'\n';
            detectionLatency << calcDetectLatencies(rippleBoundStart, rippleBoundEnd, detectionTimeIndexes, bootstrapStartIntervals[xa], bootstrapStartIntervals[xa]+minIndexCount, blockLength) <<'\n';
            relativeDetectionLatency << calcRelativeDetectLatencies(rippleBoundStart, rippleBoundEnd, detectionTimeIndexes, bootstrapStartIntervals[xa], bootstrapStartIntervals[xa]+minIndexCount, blockLength) <<'\n';
            tpRate.flush(); fpRate.flush(); detectionLatency.flush(); fpPercentage.flush();
        }
    }
    return 0;
}
