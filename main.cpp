#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <stdio.h>
#include <cstring>
#include <iomanip>
#include <pthread.h>
#include <time.h>
#include <sys/time.h>

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}



using namespace std;

void multc (long double* a, long double* b, long double* r);
void addc (long double* a, long double* b, long double* r);
void subc (long double* a, long double* b, long double* r);
long int check_iterations (long double* c, long int depth);
unsigned int* rcolor (long unsigned int j);

void *distributed_calc (void * arg);

long int** dot;
long int depth;
long unsigned int resolution;
long double sxc;
long double syc;
long unsigned int * thread_distribution;
long double * thread_spoints;
long double hsize;
long double vsize;
int threadnum;
pthread_mutex_t mut;
bool * tbcalc;
int partsnum;

int main(int argc, const char * argv[]) {
    string input = "u";
    string floc;
    bool oneshot = 0;
    unsigned int* rcol = new unsigned int [3];
    depth = 100;
    hsize = 4;
    vsize = 4;
    resolution = 500;
    sxc = -0;
    syc = -0;
    threadnum = 8;

    for (int i=0; i<argc; i++) {
        if (strcmp(argv[i], "-sxc") == 0) {
            sxc = atof(argv[i+1]);
        }
        if (strcmp(argv[i], "-syc") == 0) {
            syc = atof(argv[i+1]);
        }
        if (strcmp(argv[i], "-res") == 0) {
            resolution = atoi(argv[i+1]);
        }
        if (strcmp(argv[i], "-depth") == 0) {
            depth = atoi(argv[i+1]);
        }
        if (strcmp(argv[i], "-size") == 0) {
            vsize = atof(argv[i+1]);
            hsize = vsize;
        }
        if (strcmp(argv[i], "-fnum") == 0) {
            floc = "/home/nikita/fractals/" + to_string(atoi(argv[i+1])) + ".ppm";
        }
        if (strcmp(argv[i], "-oneshot") == 0) {
            oneshot = 1;
        }
        if (strcmp(argv[i], "-threadnum") == 0) {
            threadnum = atoi(argv[i+1]);
        }
    }

    pthread_t thread[threadnum];
    partsnum = 80;


    thread_distribution = new long unsigned int [partsnum];
    thread_spoints = new long double [partsnum];
    tbcalc = new bool [partsnum];

    input = "calc";
    while (input != "exit") {
        while (input != "calc") {
        for (int i =0; i<partsnum; i++) {
            tbcalc[i] = 0;
        }
        cin >> input;
            if(input == "a"){
                sxc = sxc - hsize/10;
            }
            if(input == "d"){
                sxc = sxc + hsize/10;
            }
            if(input == "s"){
                syc = syc + vsize/10;
            }
            if(input == "w"){
                syc = syc - vsize/10;
            }
            if(input == "+"){
                vsize = vsize*2/3;
                hsize = hsize*2/3;
            }
            if(input == "-"){
                vsize = vsize/(2*1.0/3);
                hsize = hsize/(2*1.0/3);
            }
            if(input == "hr"){
                resolution = 1000;
            }
            if(input == "lr"){
                resolution = 500;
            }
            if(input == "ehr"){
                resolution = 3000;
            }
            if(input == "ldepth"){
                depth = 100;
            }
            if(input == "mdepth"){
                depth = 150;
            }
            if(input == "hdepth"){
                depth = 200;
            }
            if(input == "ehdepth"){
                depth = 300;
            }
            if(input == "d++"){
                depth+=50;
            }
            if(input == "d--"){
                depth-=50;
            }
            if(input == "r++"){
                resolution+=500;
            }
            if(input == "r--"){
                resolution-=500;
            }
            if(input == "exit"){
                return 0;
            }
            if(input == "r+"){
                resolution+=10;
            }
            if(input == "r-"){
                resolution-=10;
            }
            if(input == "setr"){
                cin >> resolution;
            }
            if(input == "setd"){
                cin >> depth;
            }
            if(input == "pp"){
                cout << endl << "///////////////////////";
                cout << endl <<"resolution:" << endl;
                cout << resolution << endl << endl;
                cout << "depth:" << endl;
                cout << depth << endl << endl;
                cout << "x coord:" << endl;
                cout << setprecision(40) << sxc << endl << endl;
                cout << "y coord:" << endl;
                cout << setprecision(40) << syc << endl << endl;
                cout << "size:" << endl;
                cout << setprecision(40) << vsize << endl << endl;
                cout << "///////////////////////" << endl;
            }

        }

        dot = new long int * [resolution];
        for(long unsigned int i=0; i< resolution; i++){
            dot[i] = new long int [resolution];
        }
        for (int i=0; i<partsnum; i++) {
            thread_distribution[i] = int(((i+1)*resolution)/(partsnum));
            thread_spoints[i] = syc - vsize/2 + (i*vsize)/partsnum;

        }
        for (int i =0; i<threadnum; i++) {
            pthread_create(&thread[i], NULL, distributed_calc, (void *) i);
        }
        //time_t t = time(NULL);
        double t = get_wall_time();
        cout << "waiting for other threads" << endl;
        for (int i=0; i<threadnum; i++) {
            pthread_join(thread[i],NULL);
        }

        t = get_wall_time() - t;
        cout << "All threads finished in " << t << " seconds" << endl;
        cout << floc;


        FILE *fp = fopen(floc.c_str(), "wb");
        (void) fprintf(fp, "P6\n%d %d\n255\n", resolution, resolution);

        for (long unsigned int i = 0; i < resolution; i++)
        {
            for (long unsigned int j = 0; j < resolution; j++)
            {
                static unsigned char color[3];
                rcol = rcolor(dot[i][j]);
                color[0] = rcol[0];
                color[1] = rcol[1];
                color[2] = rcol[2];
                //255*cos((((dot[i][j]*1.0)/50)*((dot[i][j]*1.0)/50))*(3.14/2))
                (void) fwrite(color, 1, 3, fp);
            }
            if (i == resolution - 1) {
                cout << endl << "Finished" << endl;
            }
        }
        (void) fclose(fp);


        for( long int i=0; i< resolution; i++) {
            delete [] dot[i];
        }
        delete [] dot;
        input = "u";
        if (oneshot){
        break;
        }
    }
}

void *distributed_calc (void * arg) {
    long double* coord;
    long unsigned int spoint;
    int tnum = 0;
    choose_sector:
    tnum = 0;
    pthread_mutex_lock(&mut);
    while (tnum < partsnum ) {
        if (tnum == partsnum-1 && tbcalc[tnum] == 1) {
            goto finish;
        }
        if (tbcalc[tnum] == 0) {
            tbcalc[tnum] = 1;
            break;
        }
        tnum++;
    }
    pthread_mutex_unlock(&mut);

    coord = new long double [2];
    coord[1] = thread_spoints[tnum];
    coord[0] = sxc-vsize/2;

    if (tnum == 0) {
        spoint = 0;
    } else {
        spoint = thread_distribution[tnum-1];
    }



    for (long unsigned int i=spoint; i<thread_distribution[tnum]; i++){

        for (long unsigned int j=0; j<resolution; j++){


            dot[i][j] = check_iterations(coord, depth);


            coord[0] += vsize/resolution;
            if(i==thread_distribution[tnum]-1 && j == resolution-1) {
                pthread_mutex_lock(&mut);
                cout << "sector №" << tnum << " of " << partsnum << " finished"<<endl;
                pthread_mutex_unlock(&mut);
            }
        }
        coord[1] += hsize/resolution;
        coord[0] = sxc-vsize/2;


    }
    delete [] coord;
    goto choose_sector;
    finish:
    pthread_mutex_unlock(&mut);
}

void func (long double* a, long double* r, long double* c) {
    long double* tmp = new long double [2];
    long double* tmp2 = new long double [2];
    multc(a,a,tmp);
    addc(tmp,c,tmp2);
    r[0] = tmp2[0];
    r[1] = tmp2[1];
    delete [] tmp;
    delete [] tmp2;
    return;
}

void multc (long double* a, long double* b, long double* r) {
    r[0] = a[0]*b[0] - a[1]*b[1];
    r[1] = a[0]*b[1] + a[1]*b[0];
    return;
}

void addc (long double* a, long double* b, long double* r) {
    r[0] = a[0] + b[0];
    r[1] = a[1] + b[1];
    return;
}

void subc (long double* a, long double* b, long double* r) {
    r[0] = a[0] - b[0];
    r[1] = a[1] - b[1];
    return;
}

long int check_iterations (long double* c, long int depth) {
    long int niter = 0;
    long double* a = new long double [2];
    long double* r = new long double [2];
    a[0] = 0;
    a[1] = 0;
    while ( a[0]*a[0] + a[1]*a[1] <=9 ) {
        niter++;
        func(a,r,c);
        a[0] = r[0];
        a[1] = r[1];
        if (niter == depth){
            niter = -1;
            break;
        }
    }
    delete [] a;
    delete [] r;
    return niter;
}

unsigned int* rcolor (long unsigned int j) {
    int cnum = 3;
    long unsigned int res = depth;
    j = j % res;
    unsigned int** r = new unsigned int* [cnum];

    for (int i=0; i<cnum; i++){
        r[i] = new unsigned int [3];
    }
    r[0][0] = 255;
    r[0][1] = 255;
    r[0][2] = 0;

    r[1][0] = 0;
    r[1][1] = 0;
    r[1][2] = 255;

    r[2][0] = 0;
    r[2][1] = 0;
    r[2][2] = 0;


    long unsigned int* distrib = new long unsigned int [cnum];
    for (int i=1; i<cnum-1; i++) {
        distrib[i] = int(i*res/(cnum-1));
    }
    distrib[0] = 0;
    distrib[cnum-1] = res;

    unsigned int* rtrn = new unsigned int [3];
    int pl = 0;

    for (int i=0; i<cnum; i++) {
        if(i!=cnum-1) {
            if(j >= distrib[i] && j <= distrib[i+1]) {
                pl = i;
            }
        }
    }
    for(int i=0; i<3; i++) {
        rtrn[i] = int( (r[pl][i]*(1-(j-distrib[pl]))/(distrib[pl+1]-distrib[pl])) + (r[pl+1][i]*((distrib[pl+1]-j))/(distrib[pl+1]-distrib[pl])) );
    }
    return rtrn;
}
