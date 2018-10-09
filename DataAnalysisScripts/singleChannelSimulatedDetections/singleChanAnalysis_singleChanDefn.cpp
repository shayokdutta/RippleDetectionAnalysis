#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <sstream>
#include <pthread.h>
#include <unistd.h>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

#define THRESHOLDSTART 6.25
#define THRESHOLDEND 10.25
#define THRESHOLDSET false
#define THRESHOLDT2 362

//perform simulated real-time ripple detection on this file
#define DATAINT2FILENAME "/home/shayok/Documents/Code/RippleDetectionAnalysis/Cavaradossi/paperData/singleChanAnalysis/smoothed_envelope_simulated.out"

#define BOOTSTRAPS 1000

//canonical start and end times
#define RIPPLESTARTBOUND "/home/shayok/Documents/Code/RippleDetectionAnalysis/Cavaradossi/paperData/singleChanAnalysis/rippleBoundsStart3SD.out"
#define RIPPLEENDBOUND "/home/shayok/Documents/Code/RippleDetectionAnalysis/Cavaradossi/paperData/singleChanAnalysis/rippleBoundsEnd3SD.out"

//output files *threshold extentions added within analysis code
#define SIMDETECTIONFILENAME "/home/shayok/Documents/Code/RippleDetectionAnalysis/Cavaradossi/paperData/offlineAnalysis/singleChan/singleChanDefn3SD/simDetectionsSingleChan"
#define TPRATEFILENAME "/home/shayok/Documents/Code/RippleDetectionAnalysis/Cavaradossi/paperData/offlineAnalysis/singleChan/singleChanDefn3SD/tpRate"
#define FPRATEFILENAME "/home/shayok/Documents/Code/RippleDetectionAnalysis/Cavaradossi/paperData/offlineAnalysis/singleChan/singleChanDefn3SD/fpRate"
#define FPPERCENTAGEFILENAME "/home/shayok/Documents/Code/RippleDetectionAnalysis/Cavaradossi/paperData/offlineAnalysis/singleChan/singleChanDefn3SD/fpPercent"
#define DETECTIONLATENCYFILENAME "/home/shayok/Documents/Code/RippleDetectionAnalysis/Cavaradossi/paperData/offlineAnalysis/singleChan/singleChanDefn3SD/detectionLatency"
#define RELATIVEDETECTIONLATENCYFILENAME "/home/shayok/Documents/Code/RippleDetectionAnalysis/Cavaradossi/paperData/offlineAnalysis/singleChan/singleChanDefn3SD/relativeDetectionLatency"

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

//file handlers struct for further multithreading if we want to debug faster...
struct filehandlers{
    std::ofstream tpRate;
    std::ofstream fpRate;
    std::ofstream detectionLatency;
    std::ofstream relativeDetectionLatency;

    std::ifstream startBounds;
    std::ifstream endBounds;
    std::ifstream dataInT2;
};


void* real_work_thread(void *arg)
{
    double x = *((double*) arg); //stores threshold
    std::cout<<"Single Chan Detection Threshold " <<x<<" thread running!"<<'\n';

    //Read input files
    std::string line;
    struct filehandlers fileHandlers;
    fileHandlers.startBounds.open(RIPPLESTARTBOUND);
    std::vector <int> rippleBoundStart;
    if(fileHandlers.startBounds.is_open()){
        while(std::getline(fileHandlers.startBounds,line)){
            rippleBoundStart.push_back(std::stod(line));
        }
    }
    else{
        std::cout<< "Error opening file ripple bounds start"<<std::endl;
        pthread_exit(0);
    }
    fileHandlers.startBounds.close();

    fileHandlers.endBounds.open(RIPPLEENDBOUND);
    std::vector <int> rippleBoundEnd;
    if(fileHandlers.endBounds.is_open()){
        while(std::getline(fileHandlers.endBounds,line)){
            rippleBoundEnd.push_back(std::stod(line));
        }
    }
    else{
        std::cout<< "Error opening file ripple bounds end"<<std::endl;
        pthread_exit(0);
    }
    fileHandlers.endBounds.close();

    fileHandlers.dataInT2.open(DATAINT2FILENAME);
    std::vector <double> smoothed_envelopeT2;
    if(fileHandlers.dataInT2.is_open()){
        while(std::getline(fileHandlers.dataInT2,line)){
            smoothed_envelopeT2.push_back(std::stod(line));
        }
    }
    else{
        std::cout<<"Error opening file T2"<<std::endl;
        pthread_exit(0);
    }
    fileHandlers.dataInT2.close();

    int blockLength = 600; //200ms block after detection or 600 indexes at sampling rate 3kHz
    //open output detection file
    std::ostringstream strs;
    strs << x*100;
    std::ofstream myfile;
    std::string fileNameee = SIMDETECTIONFILENAME;
    fileNameee += strs.str() + ".out";
    myfile.open(fileNameee, std::ofstream::out | std::ofstream::trunc);

    //detect and store!
    std::vector <int> detectionTimeIndexes;
    //set detection threshold
    double thresholdT2;
    if(THRESHOLDSET)
            thresholdT2 = THRESHOLDT2;
        else
            thresholdT2 = calcMean(smoothed_envelopeT2)+(x*calcSTD(smoothed_envelopeT2));

    

    //Hunt for ripples 
    unsigned int i=0;
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

    //track which detections detect canonical ripples
    std::vector<std::string> rippleDetected;

    //all detections initially marked as not been detected
    for(unsigned int xx = 0; xx<rippleBoundStart.size(); ++xx){
        rippleDetected.push_back("F");
    }
    unsigned int yyy = 0;
    unsigned int iiiii = 0;
    while(iiiii<detectionTimeIndexes.size() && yyy < rippleBoundStart.size()){
        //if pre detect!
        if(detectionTimeIndexes[iiiii]<rippleBoundStart[yyy]){
            //don't double penalize for both false detection and missed detection
            if(detectionTimeIndexes[iiiii]+blockLength > rippleBoundStart[yyy]){
                /**
                 * It's worth noting here that this accounts for the case that
                 * there are 3 ripples consecutively and the first one is the 
                 * one that is detected. The following one is missed because it 
                 * is within 200 ms of the first one (so practically it would've
                 * been disrupted as well). The second one would be missed and
                 * the third one would be detected but it may be a late detection
                 * thus penalizing the detection latency quantification. Skipping
                 * events like this (they do occur occassionally) provides a more
                 * representative view of boths system and algorithm performance.
                 */
               rippleDetected[yyy] = "skip"; ++yyy;
            }
            ++iiiii;
        }
        //if ripple detected!
        else if(detectionTimeIndexes[iiiii] >= rippleBoundStart[yyy] 
            && detectionTimeIndexes[iiiii] <= rippleBoundEnd[yyy]){
                rippleDetected[yyy]=std::to_string(iiiii);
                //make sure we don't penalize detection that would be disrupted
                if(yyy+1 < rippleDetected.size()){
                    if(rippleBoundStart[yyy+1]-rippleBoundEnd[yyy]<blockLength){
                        rippleDetected[yyy+1]="skip";
                        ++yyy;
                    }
                }
                else
                    break;
                ++iiiii;++yyy;
        }
        //if ripple was missed!
        else if(detectionTimeIndexes[iiiii]>rippleBoundEnd[yyy])
        {
            //std::cout<<"MISSED DETECTION RB INDEX: " <<yyy<<'\n';
            ++yyy;
        }
    }

    //open detection metric quantification files
    
    fileNameee = TPRATEFILENAME;
    fileNameee += strs.str() + ".out";
    fileHandlers.tpRate.open(fileNameee, std::ofstream::out | std::ofstream::trunc);

    fileNameee = FPRATEFILENAME;
    fileNameee += strs.str() + ".out";
    fileHandlers.fpRate.open(fileNameee, std::ofstream::out | std::ofstream::trunc);

    fileNameee =DETECTIONLATENCYFILENAME;
    fileNameee += strs.str() + ".out";
    fileHandlers.detectionLatency.open(fileNameee, std::ofstream::out | std::ofstream::trunc);


    fileNameee =RELATIVEDETECTIONLATENCYFILENAME;
    fileNameee += strs.str() + ".out";
    fileHandlers.relativeDetectionLatency.open(fileNameee, std::ofstream::out | std::ofstream::trunc);


    //Generate 1000 full dataset sample metric quantifications
    /**
     * NOTE: THIS IS NOT GOING TO BE OPTIMIZED I'M SIMPLY MODIFYING THE 20 MIN
     * CHUNK VERSION TO FULL DATASET AND LEAVING THIS RUNNING WHILE I DO OTHER
     * STUFF (e.g., eat lunch, get my microdrive ready)
     */
    //chunk data
    std::vector<int> bootstrapsampleStartTimes;
    unsigned int myCounter = 0;
    while(myCounter < smoothed_envelopeT2.size()){
        bootstrapsampleStartTimes.push_back(myCounter);
        myCounter+=45;//take every 45th sample so every 15ms...this is binning
    }
    myCounter = 0;
    auto t1 = Clock::now();
    while(myCounter < BOOTSTRAPS){
        ++myCounter;
        if(myCounter %10 == 0){
            std::cout<<"Single Chan Detection Thread: " << x << ": Last 10 iterations time: " <<
            std::chrono::duration_cast<std::chrono::seconds>(Clock::now() - t1).count()
            <<" s" << ". Percent Complete : "<<(double)myCounter*100/BOOTSTRAPS<<"%"<<std::endl;
            t1 = Clock::now();
        }
        //Metric quantification algorithm
        //instantiate detection latency metrics vectors
        std::vector<double> detectionlatency, relativedetectionlatency;
        //instantiate canonical ripple counter
        int totalRipples = 0;
        //instantiate true/false positive counters
        int trueDetections = 0;
        int falseDetections = 0;
        //instantiate true negative counter
        int totalTrueNegatives = 0;

        // auto timeSection = Clock().now();

        //Generate full dataset bootstrap samples
        for(unsigned int iii=0; iii<bootstrapsampleStartTimes.size(); ++iii){
            //generate random integer from 0 to size of data
            int bootstrapsample = bootstrapsampleStartTimes[std::rand() % (bootstrapsampleStartTimes.size())];
            
            /** Determine if sample is a canonical ripple
             * tail of 15 ms chunk has ripple bound start
             * 15 ms chunk is within ripple bound start and end
             * beginning of 15 ms chunk has ripple bound end
             *  this isn't a conditional b/c we're only looking at start times
             *  it's already accounted for by the  falling within the start and end
             * store associated ripple bound start and end if ripple
            */
            //Determine if sample is a canonical ripple
            int rippleBoundSample = -1;
            //loop through ripple bounds
            for(unsigned int xx=0; xx<rippleBoundStart.size(); ++xx){
                //if 15 ms chunk is within ripple start and end
                if(bootstrapsample > rippleBoundStart[xx] && bootstrapsample < rippleBoundEnd[xx]){
                    rippleBoundSample = xx;
                    break;
                }
                //if tail of 15 ms sample has ripple bound start
                else if(bootstrapsample+45 > rippleBoundStart[xx] && bootstrapsample+45 < rippleBoundEnd[xx]){
                    rippleBoundSample = xx;
                    break;
                }
            }
            //if sample is canonical ripple and we don't want to skip it because of block length
            if(rippleBoundSample != -1 && 
                rippleDetected[rippleBoundSample] != "skip"){
                //increment canonical ripple counter
                ++totalRipples;
                //Determine if canonical ripple was detected
                if(rippleDetected[rippleBoundSample].compare("F") != 0 &&
                    rippleDetected[rippleBoundSample].compare("skip") != 0){
                    int xx = 0;
                    std::stringstream streeem(rippleDetected[rippleBoundSample]);
                    streeem >> xx;
                    ++trueDetections;
                    //quantify detection latency in ms and store in vector (ms)
                    double detectionlat = (detectionTimeIndexes[xx]-rippleBoundStart[rippleBoundSample])/3;
                    detectionlatency.push_back(detectionlat);
                    //quantify relative detection latency and store in vector
                    relativedetectionlatency.push_back(detectionlat/((rippleBoundEnd[rippleBoundSample]-rippleBoundStart[rippleBoundSample])/3));
                }
            }
            //else determine if detection is present for false detection
            else if(rippleBoundSample == -1){
                ++totalTrueNegatives;
                //loop through all simulated detection time indices
                for(unsigned int xx=0; xx<detectionTimeIndexes.size(); ++xx){
                    //if detection is present within random sample bound
                    if(detectionTimeIndexes[xx] > bootstrapsample 
                        && detectionTimeIndexes[xx] < bootstrapsample + 45){
                            //increment false positive counter
                            ++falseDetections;
                            break;
                        }
                }
            }
        }
        /** Calculate metrics
         * tprate = tpcount/total canonical ripples
         * fprate = fpcount/20 min
         * average detection latencies
        */
        //write above metrics to appropriate files.
        fileHandlers.tpRate << (double)trueDetections/totalRipples << '\n';
        fileHandlers.fpRate << (double)falseDetections/(totalTrueNegatives*15/60000)<< '\n';
        fileHandlers.detectionLatency << calcMean(detectionlatency)<< '\n';
        fileHandlers.relativeDetectionLatency << calcMean(relativedetectionlatency)<< '\n';
        fileHandlers.tpRate.flush();fileHandlers.fpRate.flush();fileHandlers.detectionLatency.flush();fileHandlers.relativeDetectionLatency.flush();
        usleep(10);
    }
    std::cout<<"Threshold: " <<x<<" thread complete!!"<<'\n';
    fileHandlers.tpRate.close();fileHandlers.fpRate.close();fileHandlers.detectionLatency.close();fileHandlers.relativeDetectionLatency.close();
    pthread_exit(0);
}

int main(int argc, char *argv[])
{
    std::vector<pthread_t> tids;
    for(double x = THRESHOLDSTART; x<THRESHOLDEND; x+=0.25){ //loop through a bunch of thresholds, perform ripple detections, and evaluate
        pthread_t tid;
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_create(&tid,&attr,real_work_thread,&x);
        tids.push_back(tid);
        sleep(1);
    }
    for(unsigned int x = 0; x<tids.size(); ++x){
        int rc = pthread_join(tids[x],NULL);
        if (rc) { fprintf(stderr, "failed to join thread #%lf - %s\n",
                                (long)(x*0.25)+1.5, strerror(rc));
               exit(EXIT_FAILURE);}
    }
    return 0;
}
