#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <queue>
#include <fstream>

// ----------------- SIMULATION PARAMETERS -----------------
// SELECT ONE OF THE FOLLOWING DAAS (NEFDA or MODIFIED_NEFDA) TO RUN THE SIMULATION
// #define NEFDA
#define MODIFIED_NEFDA

// SELECT ONE OF THE FOLLOWING HASH RATE FUNCTIONS (NON_ADAPTIVE OR ADAPTIVE) TO RUN THE SIMULATION
#define NON_ADAPTIVE
// #define ADAPTIVE
// ----------------- SIMULATION PARAMETERS -----------------

#define NETWORK_RESOLUTION 1.0

#define DAY (86400.0)
#define FORTNIGHT (14.0*DAY)
#define YEAR (365.25*DAY)
#define INITIAL_ABS_DIFFICULTY (4294967296.0) // 2^32
// #define INITIAL_ABS_DIFFICULTY (pow(2,69)) // ~ 2^69
#define INITIAL_HASH_RATE (10000000.0)
// #define INITIAL_HASH_RATE (pow(10,18))


// standard NEFDA algorithm
#ifdef NEFDA
#define SHADOW_OF_PAST_TIMESCALE (FORTNIGHT)
#endif

// modified NEFDA algorithm
#ifdef MODIFIED_NEFDA
#define SHADOW_OF_PAST_TIMESCALE (FORTNIGHT/0.5)
#define SHADOW_OF_FUTURE_TIMESCALE (FORTNIGHT/0.5)
#endif

#define IDEAL_INTERBLOCK_TIME (600.0)
#define GENESIS_TIME (2009.0*YEAR)
#define FINISH_TIME (2024.0*YEAR)
#define NUM_OF_LAST_BLOCKS 6

using namespace std;

typedef struct {
  double t;
  double target;
  double faded_total_work;
} BlockTemplate;


#ifdef NON_ADAPTIVE
double non_adaptive_hash_rate_function(double t) {
    // constant hash function
    // return INITIAL_HASH_RATE;

    // step hash function
    // return (t < 2010.5*YEAR) ? pow(10,18) : pow(10,18)*30;

    // step hash function
    // return (t < 2010.5*YEAR) ? 10000000 : 10000000*30;

    // negative step hash function
    // return (t < 2010.5*YEAR) ? 10000000 : 10000000/30;

    // cos hash function
    // return (10000000 - 5000000 * sin(2*M_PI*(t-GENESIS_TIME)/(4 * YEAR)));

    // rectangular hash function
    // if (t < 2010.5*YEAR) return 10000000;
    // else if(t < 2016.5*YEAR) return 10000000*30;
    // else return 10000000;

    // staircase hash function
    // if (t < 2010.5*YEAR) return 10000000;
    // else if(t < 2013.5*YEAR) return 10000000*10;
    // else return 10000000 * 100;

    // linearly increasing hash function
    // return INITIAL_HASH_RATE + (INITIAL_HASH_RATE/YEAR) * (t-GENESIS_TIME);

    // quadratically increasing hash function
    // return INITIAL_HASH_RATE + INITIAL_HASH_RATE/pow(YEAR,2) * pow((t-GENESIS_TIME),2);

    // exponentially increasing hash function
    // return INITIAL_HASH_RATE * exp(log(2) * (t-GENESIS_TIME)/YEAR);

    // superexponentially increasing hash function
    // return INITIAL_HASH_RATE * exp(log(2) * pow(((t-GENESIS_TIME)/YEAR),2));
    
    // sin and linearly increasing hash function
    return INITIAL_HASH_RATE + INITIAL_HASH_RATE/YEAR * (t-GENESIS_TIME) + INITIAL_HASH_RATE * sin(2*M_PI*(t-GENESIS_TIME)/(YEAR));
}
#endif


#ifdef ADAPTIVE
double adaptive_hash_rate_function(double avg_difficulty){
    double hb = INITIAL_HASH_RATE / 3;
    double hv = 4 * hb;
    double hg = 4 * hb;
    double h = hb;
    double epsilon = 0.03;

    if(avg_difficulty <= 1 - epsilon){
        h += hg + hv;
    }
    else if(avg_difficulty <= 1 - epsilon/3){
        h += hg + hv / (1 + exp((-1)*(1-avg_difficulty)/(epsilon*epsilon)));
    }
    else if(avg_difficulty <= 1 + epsilon){
        h += hv / (1 + exp((-1)*(1-avg_difficulty)/(epsilon*epsilon)));
    }
    return h;
}
#endif


int main() {

    FILE *fp;
    #ifdef NEFDA
        fp = fopen("data/NEFDA_Output.csv", "w");
    #endif
    #ifdef MODIFIED_NEFDA
        fp = fopen("data/MODIFIED_NEFDA_Output.csv", "w");
    #endif

    srand48(time(NULL));

    #ifdef ADAPTIVE
    // queue for avg difficulty/avg target
    queue < double > q;
    double running_sum = NUM_OF_LAST_BLOCKS;
    double running_avg = 1;
    for (int i = 0; i < NUM_OF_LAST_BLOCKS; i++){
        q.push(1.0);
    }
    #endif
    
    BlockTemplate *blocks = (BlockTemplate *) calloc(1000000, sizeof(BlockTemplate));
    int n=0;
    double t = GENESIS_TIME;  
    BlockTemplate candidate = {t, 1.0/INITIAL_ABS_DIFFICULTY, INITIAL_ABS_DIFFICULTY*(SHADOW_OF_PAST_TIMESCALE/IDEAL_INTERBLOCK_TIME)};
    blocks[n] = candidate; // genesis block

    fprintf(fp, "BLOCK,timestamp,t,block_time,real_hash_rate,difficulty,faded_total_work,skew\n");	 
    #ifdef NON_ADAPTIVE      
    fprintf(fp, "%d,%.0f,%.10lf,%.10lf,%.0f,%.10lf,%.10lf,%0.10lf\n", n, blocks[n].t, blocks[n].t/YEAR, -1.0, non_adaptive_hash_rate_function(t), 1.0/(blocks[n].target)/INITIAL_ABS_DIFFICULTY, blocks[n].faded_total_work, (blocks[n].t-(GENESIS_TIME+n*IDEAL_INTERBLOCK_TIME))/IDEAL_INTERBLOCK_TIME);	       
    #endif
    #ifdef ADAPTIVE
    fprintf(fp, "%d,%.0f,%.10lf,%.10lf,%.0f,%.10lf,%.10lf,%0.10lf\n", n, blocks[n].t, blocks[n].t/YEAR, -1.0, adaptive_hash_rate_function(running_avg), 1.0/(blocks[n].target)/INITIAL_ABS_DIFFICULTY, blocks[n].faded_total_work, (blocks[n].t-(GENESIS_TIME+n*IDEAL_INTERBLOCK_TIME))/IDEAL_INTERBLOCK_TIME);	       
    #endif
        // BLOCK 0                                      // block number
        // t=2009.0000000000(-1.0000000000)             // current time (block time = blocks[n].t-blocks[n-1].t) 
        // diff=1.0000000000                            // difficulty = 1 / (blocks[n].target * INITIAL_ABS_DIFFICULTY)
        // ftw=8658654068736.0000000000                 // faded total work = blocks[n].faded_total_work 
        // skew=0.0000000000                            // skew = (blocks[n].t - (GENESIS_TIME+n*IDEAL_INTERBLOCK_TIME)) / IDEAL_INTERBLOCK_TIME
                                                        // skew = (current time - (GENESIS_TIME + block number * IDEAL_INTERBLOCK_TIME)) / IDEAL_INTERBLOCK_TIME
                                                        // skew = drift of current time from the ideal time / ideal block time
                                                        // skew = how many ideal block times are we ahead or behind
                                                        // of course here skew is 0 because we are at the genesis block

    while ( (n + 1) < 1000000 && t < FINISH_TIME) {     // repeat until we reach 1000000 blocks or finish time

        do {                                            // procedure that for every t (in intervals of NETWORK_RESOLUTION) mints a candidate block, breaks if the candidate block is valid
            t += NETWORK_RESOLUTION;
            candidate.t = t;
            double estimated_faded_total_work = blocks[n].faded_total_work*exp(-(t-blocks[n].t)/SHADOW_OF_PAST_TIMESCALE);      
            double estimated_hash_rate = estimated_faded_total_work/SHADOW_OF_PAST_TIMESCALE;

            #ifdef NEFDA
                candidate.target = 1.0/(estimated_hash_rate*IDEAL_INTERBLOCK_TIME);
            #endif

            #ifdef MODIFIED_NEFDA
                double skew = t-(GENESIS_TIME+(n+1)*IDEAL_INTERBLOCK_TIME);
                candidate.target = 1.0/(estimated_hash_rate*IDEAL_INTERBLOCK_TIME)*exp(skew/SHADOW_OF_FUTURE_TIMESCALE);
            #endif
            
            candidate.faded_total_work = estimated_faded_total_work + 1.0/(candidate.target);
        #ifdef NON_ADAPTIVE
        } while (!(drand48() < candidate.target*non_adaptive_hash_rate_function(t)*NETWORK_RESOLUTION ));      // if the candidate block is valid, break (the mul with non_adaptive_hash_rate_function(t) simulates the effect of volatile network hash rate)
        #endif
        #ifdef ADAPTIVE
        } while (!(drand48() < candidate.target*adaptive_hash_rate_function(running_avg)*NETWORK_RESOLUTION ));      // if the candidate block is valid, break (the mul with adaptive_hash_rate_function(t) simulates the effect of volatile network hash rate)
        #endif

        n++;
        blocks[n] = candidate;                              // add the new valid block to the blocks array  
        
        #ifdef NON_ADAPTIVE
        fprintf(fp, "%d,%.0f,%.10lf,%.10lf,%.0f,%.10lf,%.10lf,%0.10lf\n", n, blocks[n].t, blocks[n].t/YEAR, blocks[n].t-blocks[n-1].t, non_adaptive_hash_rate_function(t), 1.0/(blocks[n].target)/INITIAL_ABS_DIFFICULTY, blocks[n].faded_total_work, (blocks[n].t-(GENESIS_TIME+n*IDEAL_INTERBLOCK_TIME))/FORTNIGHT);	
        #endif
        #ifdef ADAPTIVE
        fprintf(fp, "%d,%.0f,%.10lf,%.10lf,%.0f,%.10lf,%.10lf,%0.10lf\n", n, blocks[n].t, blocks[n].t/YEAR, blocks[n].t-blocks[n-1].t, adaptive_hash_rate_function(running_avg), 1.0/(blocks[n].target)/INITIAL_ABS_DIFFICULTY, blocks[n].faded_total_work, (blocks[n].t-(GENESIS_TIME+n*IDEAL_INTERBLOCK_TIME))/FORTNIGHT);	
        #endif
            // BLOCK 796672                                 // block number
            // t=2023.9999586154(905.0000000000)            // current time (block time in NETWORK_RESOLUTION i.e. seconds = blocks[n].t-blocks[n-1].t)
            // diff=42.0106884009                           // difficulty = 1 / (blocks[n].target * INITIAL_ABS_DIFFICULTY)
            // ftw=363936452585549.8750000000               // faded total work = blocks[n].faded_total_work
            // skew=-3.8363971561                           // skew = (blocks[n].t - (GENESIS_TIME+n*IDEAL_INTERBLOCK_TIME)) / FORTNIGHT
                                                            // skew = (current time - (GENESIS_TIME + block number * IDEAL_INTERBLOCK_TIME)) / FORTNIGHT
                                                            // skew = drift of current time from the ideal time / FORTNIGHT
                                                            // skew = how many FORTNIGHTS are we ahead or behind
                                                            // here skew is -3.8363971561 because we are 3.8363971561 FORTNIGHTS ahead of the ideal time

        #ifdef ADAPTIVE
        running_sum += 1.0/(blocks[n].target)/INITIAL_ABS_DIFFICULTY;
        q.push(1.0/(blocks[n].target)/INITIAL_ABS_DIFFICULTY);
        running_sum = running_sum - q.front();
        q.pop();
        running_avg = running_sum / NUM_OF_LAST_BLOCKS;
        #endif
    }

    free(blocks);
    fclose(fp);
    return 0;
}

