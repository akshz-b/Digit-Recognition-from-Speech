// 234101006_digit_recognition.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <float.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <vector>

#include <algorithm>
#include <sstream>

#define CONVERGE_ITERATIONS 200
#define M 32 //Number of obsevation symbols per state
#define N 5 //Number of states
#define P 12
#define LIMIT 5000
#define CB_SIZE 32
#define PI 3.142857142857
#define FRAME_SIZE 320
#define FRAMES 100

const int WIN_SIZE  = (FRAME_SIZE * FRAMES);

int T; //Time sequence length
const int MAX_T = 150; // max time sequence length
using namespace std;

//Global variables
long double dcShift, nFactor, mx, silenceEnergy;
long double const threshold = 1e-30;   //Min threshold to be assigned to zero values in matrix B.
long int sSize = 0, sampSize = 0, enSize = 0;
long double max_pobs_model = 0;
int test_ans = 0, fake = 0;

//Globally defined arrays
int O[MAX_T+1];	//Observation sequence
int Q[MAX_T+1];	//state sequence.
long double pstar = 0, prev_p_star = -1;
long double Alpha[MAX_T+1][N+1];
long double Beta[MAX_T+1][N+1];
long double Gamma[MAX_T+1][N+1];
long double Delta[MAX_T+1][N+1];
int Psi[MAX_T+1][N+1];
long double Xi[MAX_T+1][N+1][N+1];

long double codeBook[CB_SIZE][P];

long int sample[1000000];
long double steadySamp[WIN_SIZE];
long double energy[100000];
long double Ai[P+1], Ri[P+1], Ci[P+1];
//tokhura weights
double tokhuraWeights[]={1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

//Model parameters A, B and Pi
long double A[N+1][N+1] = {0};
long double B[N+1][M+1] = {0};
long double Pi[N+1] = {0};

long double Abar[N+1][N+1] = {0};
long double Bbar[N+1][M+1] = {0};
long double Pibar[N+1] = {0};

long double a_average[N+1][N+1] = {0};
long double b_average[N+1][M+1] = {0};
long double pi_average[N+1] = {0};

//files
char* A_file = "a_i_j.txt";
char* B_file = "b_i_j.txt";
char* PI_file = "pi.txt";

int cnt = 1, train = 0;
long double P_O_given_Model = 0;
ofstream uni, dump;
FILE *common_dump;
int environment_known = 0, is_live_testing = 0;

#include"speech_features.h"

//Calculation of alpha variable to find the solution of problem number 1.
void forward_procedure(){
	int i , j , t;
	long double sum ;
	int index = O[1];
	P_O_given_Model = 0;

	for(i=1;i<=N;i++){
		Alpha[1][i] = Pi[i]*B[i][index];
	}

	for (t = 1; t < T; t++){
		for (j = 1; j <= N; j++){
			sum = 0;
			for (i = 1; i <= N; i++){
				sum += Alpha[t][i] * A[i][j];
			}
			Alpha[t + 1][j] = sum * B[j][O[t + 1]];
		}
	}
	for(i=1;i<=N;i++){
		P_O_given_Model = P_O_given_Model + Alpha[T][i];
	}
}

//Calculation of alpha variable to find the solution of problem number 1.
void forward_procedure(int iteration, FILE *fp = NULL){
	int i , j , t;
	long double sum ;
	int index = O[1];
	P_O_given_Model = 0;

	for(i=1;i<=N;i++){
		Alpha[1][i] = Pi[i]*B[i][index];
	}

	for (t = 1; t < T; t++){
		for (j = 1; j <= N; j++){
			sum = 0;
			for (i = 1; i <= N; i++){
				sum += Alpha[t][i] * A[i][j];
			}
			Alpha[t + 1][j] = sum * B[j][O[t + 1]];
		}
	}
	for(i=1;i<=N;i++){
		P_O_given_Model = P_O_given_Model + Alpha[T][i];
	}
	//finding where the model is matching
	if(P_O_given_Model > max_pobs_model){
		max_pobs_model = P_O_given_Model;
	 	test_ans = iteration;
	}

	//cout << "Digit:"<<iteration<<"\tP(obs/model) : " << P_O_given_Model <<endl;
	if(fp != NULL){
		fprintf(fp, "---> Digit %d ----- P(Obs/Model) : %g\n", iteration, P_O_given_Model);
		
	}
}

//function for testing with iteration as argument
void solution_to_prob1(int iteration, FILE *fp = NULL){
	if(fp == NULL)
		forward_procedure(iteration);
	else
		forward_procedure(iteration, fp);
}

//Calculation of Beta values
void backward_procedure(){
	int i , j , t;
	long double sum;
	int index = 0;
	for(i=1;i<=N;i++){
		Beta[T][i] = 1.0;
	}
	for(t=T-1;t>=1;t--){
		index = O[t+1];
		for(i=1;i<=N;i++){
			sum = 0;
			for(j=1;j<=N;j++){
				sum = sum + B[j][index]*A[i][j]*Beta[t+1][j];
			}
			Beta[t][i]=sum;
		}
	}
}

//Calculation gamma values
void get_gamma(){
	for(int t=1;t<=T;t++){
		for(int j=1;j<=N;j++){
			long double summation=0;
			for(int k=1;k<=N;k++){
				summation += Alpha[t][k] * Beta[t][k];
			}
			Gamma[t][j]=(Alpha[t][j] * Beta[t][j])/summation;
		}
	}
}

//Loading the model parameters with newly calculated values.
void update_Model_Parameters(){
    // Update Pi_bar
    for (int i = 1; i <= N; i++)
    {
        Pi[i] = Pibar[i];
    }

    // Update Aij_bar
    for (int i = 1; i <= N; i++)
    {
        for (int j = 1; j <= N; j++)
        {
            A[i][j] = Abar[i][j];
        }
    }

    // update Bjk_bar

    for (int j = 1; j <= N; j++)
    {
        for (int k = 1; k <= M; k++)
        {
            B[j][k] = Bbar[j][k];
        }
    }
}

void make_lambda_stochastic(){
	long double row_sum = 0  ,diff;
	long double row_max =0;
	int max_index = 0;

	// making A Stochastic
	// i range 1 to N
	for(int i=1; i<=N; i++){
		row_sum=0;
		row_max = 0;
		max_index =-1;		

		//j range 1 to N
		for(int j=1; j<=N; j++){
			if( Abar[i][j] > row_max )
			{	max_index = j;
				row_max = Abar[i][j];
			}

			row_sum += Abar[i][j];
		}//row sum
		if(row_sum != 1){
			Abar[i][max_index] = (row_sum > 1) ? (Abar[i][max_index] - (row_sum-1) ):(Abar[i][max_index] + (1-row_sum));
		}
		//Abar[i][max_index] -= (row_sum-1); 
		
	}

	// making B Stochastic
	// j range 1 to N
	for(int j=1; j<=N; j++){
		row_sum=0;
		row_max = 0;
		max_index =0;
		
		//k range 1 to M
		for(int k=1; k<=M; k++){

			if( Bbar[j][k] > row_max )
			{	max_index = k;	
				row_max =  Bbar[j][k];
			}

			row_sum += Bbar[j][k];

		}


		if(row_sum != 1){
			//B_bar[j][max_index] += 1 - row_sum;
			Bbar[j][max_index] = (row_sum > 1) ? (Bbar[j][max_index] - (row_sum-1) ):(Bbar[j][max_index] + (1-row_sum));
		}

		
	}

}

// solution to problem no. 3 of hmm
void re_estimation(){
	int i, j, k, t;
	long double sum1=0 , sum2 =0;
	
	// Calculate re-estimated Pi_bar
    for (int i = 1; i <= N; i++)
    {
        Pibar[i] = Gamma[1][i];
    }
	
	// Calculate Aij_bar
    // i and j range 1 to N
	for(i = 1; i<=N; i++){
		for(j = 1; j <= N; j++){
			 long double numerator = 0.0, denominator = 0.0;
            // sum over t range 1 to T
			for(t = 1; t <= T-1; t++){
				 numerator += Xi[t][i][j];
				 denominator += Gamma[t][i];
			}
			
			Abar[i][j] = (numerator / denominator);
		}
	}

	  // Calculate Bjk_bar
	for(j=1;j<=N;j++){
		int count=0;
		long double max=0;
		int ind_j=0, ind_k=0;

		for(k=1;k<=M;k++){
			sum1 =0 , sum2 =0;
			for(t=1;t<T;t++){
				sum1 = sum1 + Gamma[t][j];
				if(O[t]==k){
					sum2 = sum2 + Gamma[t][j];
				}
			}
			Bbar[j][k] = sum2/sum1;

			//finding max
			if(Bbar[j][k]>max){
				max=Bbar[j][k];
				ind_j = j;
				ind_k = k;
			}

			//updating new bij with threshold value if it is zero
			if(Bbar[j][k] == 0){
				Bbar[j][k]=threshold;
				count++;
			}
		}
		Bbar[ind_j][ind_k] = max - count*threshold;
	}
	
	make_lambda_stochastic();
	update_Model_Parameters();
}

//calculating xi
void get_xi(){

	 // t range 1 to T-1
    for (int t = 1; t <= T-1; t++)
    {
        long double denominator = 0.0;

        // Calculate denominator
        // sum over all j range 1 to N
        for (int i = 1; i <= N; i++)
        {

            // sum over all j range 1 to N
            for (int j = 1; j <= N; j++)
            {
                denominator += Alpha[t][i] * A[i][j] * B[j][O[t + 1]] * Beta[t + 1][j];
            }
        }

        // Calculate xi
        // probability of O sequence ending in state Si at time t

        // i range 1 to N
        for (int i = 1; i <= N; i++)
        {

            // j range 1 to N
            // sum over j 1 to N
            for (int j = 1; j <= N; j++)
            {   
                long double numerator = Alpha[t][i] * A[i][j] * B[j][O[t + 1]] * Beta[t + 1][j];
                Xi[t][i][j] = ((long double) numerator) / denominator;
            }
        }
    }

	
}

//Executing the Viterbi algorithm.
void viterbi(){
	// step 1 ---> Initialization

    // i range 1 to N
    for (int i = 1; i <= N; i++)
    {
        Delta[1][i] = Pi[i] * B[i][O[1]];
        Psi[1][i] = 0;
    }
    

	for (int t = 2; t <= T; t++)
    {

        // j range 1 to N
        for (int j = 1; j <= N; j++)
        {
            long double max = -1;
            int index = 0;

            // i range 1 to N
            for (int i = 1; i <= N; i++)
            {
                if (Delta[t - 1][i] * A[i][j] > max)
                {

                    // find max and change state
                    max = Delta[t - 1][i] * A[i][j];

                    // keep record of i
                    index = i;
                }
            }

            // multiply with observation symbol at j
            Delta[t][j] = max * B[j][O[t]];

            // store it in psi
            Psi[t][j] = index;
        }
    }

	

	 // step 3 ---> Termination
    
    long double max = 0;
    for(int i=1; i<=N; i++){
        if(Delta[T][i] > max) {
			max = Delta[T][i];
			Q[T] = i;
		}

        pstar = max;
    }

	// step 4 ---> State Sequence (Path) Backtracking
    // t range T-1 to 1
    for(int t = T-1; t>0; t--){
        Q[t] = Psi[t+1][Q[t+1]];
    }
}


//Save the model to a file.
void save_converged_model(FILE *fp){
	
	fprintf(fp, "---------------A Matrix----------------\n");
	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= N; j++){
			fprintf(fp,"%Le   ",A[i][j]);
		}
		fprintf(fp,"\n");
	}

	
	fprintf(fp, "---------------B Matrix----------------\n");
	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= M; j++){
			fprintf(fp, "%Le   ", B[i][j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "---------------Pi values----------------\n");
	
	for (int i = 1; i <= N; i++){
		fprintf(fp, "%Le   ", Pi[i]);
	}
	
}

//read A
bool readAMatrix(char *filename){
    fstream fin;
	fin.open(filename);

    //file does not exist
	if(!fin){
		cout<<"Couldn't open file: "<<filename<<"\n";
		return false;
	}
	long double word;

    int row = 1, col = 1;
    //until input is available
	while(fin >> word){
        col = 1;
        A[row][col++] = word;

        for(int i=2; i<=N; i++){
            fin>>word;
            A[row][col++] = word;
        }
        row++;
    }

	fin.close();
	return true;
}

//read B
bool readBMatrix(string filename){
	fstream fin;
	fin.open(filename);

    //file does not exist
	if(!fin){
		cout<<"Couldn't open file: "<<filename<<"\n";
		return false;
	}
	long double words;

    int row = 1, col = 1;

	while(fin>>words){
		col = 1;
		B[row][col++] = words;

		for(int i=1; i<M; i++){
			fin>>words;
			B[row][col++] = words;
		}
		row++;
	}
	
	fin.close();
	return true;
}

//read Pi
bool readPiMatrix(string filename){
	fstream fin;
	fin.open(filename);

    //file does not exist
	if(!fin){
		cout<<"Couldn't open file: "<<filename<<"\n";
		return false;
	}
	long double word;

    int col = 1;
    //until input is available
	while(fin >> word){
		col = 1;
        Pi[col++] = word;

        //save whole row
		for(int i=1;i<N;i++){
			fin>>word;
            Pi[col++] = word;
		}
	}

	fin.close();
	return true;
}

// Setting the model values, average model values, and bar model values to zero.
void erase_all_model(){
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			A[i][j] = 0;
			a_average[i][j] = 0;
			Abar[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			B[i][j] = 0;
			b_average[i][j] = 0;
			Bbar[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		Pi[i] = 0;
		Pibar[i] = 0;
		pi_average[i] = 0;
	}
}

//Clearing or resetting the current model.
void erase_model(){
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			A[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			B[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		Pi[i] = 0;
	}
}

// Clearing or resetting the average model.
void erase_avg_model(){
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			a_average[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			b_average[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		pi_average[i] = 0;
	}
}

//Reading the average model from a file.
void read_average_model(int digit){

	char filename[100];
	sprintf(filename, "output/avgmodels/digit_%d_A.txt", digit);
	readAMatrix(filename);

	sprintf(filename, "output/avgmodels/digit_%d_B.txt", digit);
	readBMatrix(filename);

	sprintf(filename, "output/avgmodels/digit_%d_PI.txt", digit);
	readPiMatrix(filename);
}

//Initialize the model based on the specified parameters.
void initialize_lambda(int digit, int seq, char *filename = "--"){
	char a_file[100], b_file[100], pi_file[100], obs_file[100];

	if(filename == "--"){
		readAMatrix(A_file);
		readBMatrix(B_file);
		readPiMatrix(PI_file);
	}else if(filename  == "avg"){
		read_average_model(digit);
		
	}else if(filename == "init"){
		sprintf(a_file, "validation/Digit %d/A_%d.txt", digit, digit);
		sprintf(b_file, "validation/Digit %d/B_%d.txt", digit, digit);
		sprintf(pi_file, "validation/Digit %d/pi_%d.txt", digit, digit);

		readAMatrix(a_file);
		readBMatrix(b_file);
		readPiMatrix(pi_file);
	}
}

//Incorporating the current model values into the average model.
void add_to_avg_model(){
	int i, j;
	for (i = 1; i <= N; i++){
		for (j = 1; j <= N; j++){
			a_average[i][j] += A[i][j];
		}
	}
	for (i = 1; i <= N; i++){
		pi_average[i] += Pi[i];
	}
	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= M; j++){
			b_average[i][j] += B[i][j];
		}
	}
}

// saving the final model
void save_final_model(char a_filename[], char b_filename[], char pi_filename[], int digit, int seq) {
    FILE *a_fp = fopen(a_filename, "w");
    FILE *b_fp = fopen(b_filename, "w");
    FILE *pi_fp = fopen(pi_filename, "w");

    int i, j;

    // Writing A values
    for (i = 1; i <= N; i++) {
        for (j = 1; j <= N; j++) {
            fprintf(a_fp, "%Le   ", A[i][j]);
        }
        fprintf(a_fp, "\n");
    }

    // Writing B values
    for (i = 1; i <= N; i++) {
        for (j = 1; j <= M; j++) {
            fprintf(b_fp, "%Le   ", B[i][j]);
        }
        fprintf(b_fp, "\n");
    }

    // Writing Pi values
    for (i = 1; i <= N; i++) {
        fprintf(pi_fp, "%Le   ", Pi[i]);
    }

    fclose(a_fp);
    fclose(b_fp);
    fclose(pi_fp);

    printf("Model for digit %d, utterance no. %d saved\n", digit, seq);
}

// Save the final model in the respective folder under "output/digit d/".
void save_final_model(int seq, int digit){
	char index[10];
	  char a_filename[100], b_filename[100], pi_filename[100];

	 sprintf(a_filename, "output/digit %d/model_A_%d.txt", digit, seq);
    sprintf(b_filename, "output/digit %d/model_B_%d.txt", digit, seq);
    sprintf(pi_filename, "output/digit %d/model_PI_%d.txt", digit, seq);
	 save_final_model(a_filename, b_filename, pi_filename, digit, seq);
}

// Computing the average of the average model.
void average_of_avg_model(int total_iterations){
	int i, j;
	for (i = 1; i <= N; i++){
		for (j = 1; j <= N; j++){
			a_average[i][j] /= total_iterations;

		}
	}
	for (i = 1; i <= N; i++){
		for (j = 1; j <= M; j++){
			b_average[i][j] /= total_iterations;
		}
	}
	for (i = 1; i <= N; i++){
		pi_average[i] /= total_iterations;
	}
}

// Saving the average model to a file.
void save_avg_model(int digit){
	char a_file_avg[100], b_file_avg[100], pi_file_avg[100], ind[3];

	sprintf(a_file_avg, "output/avgmodels/digit_%d_A.txt", digit);
	FILE *fp = fopen(a_file_avg, "w");
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			fprintf(fp, "%Le   ", a_average[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);


	sprintf(b_file_avg, "output/avgmodels/digit_%d_B.txt", digit);
	ofstream fout(b_file_avg);
	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			//fprintf(fp, "%Le   ", b_average[i][j]);
			fout<<b_average[i][j]<<"   ";
		}
		fout<<endl;
		//fprintf(fp, "\n");
	}
	fout.close();

	sprintf(pi_file_avg, "output/avgmodels/digit_%d_PI.txt", digit);
	fp = fopen(pi_file_avg, "w");
	for(int i=1; i<=N; i++){
		fprintf(fp, "%Le   ", pi_average[i]);
	}
	fclose(fp);
}



//Conducts training on the set of 20 files.
void training(){
	char filename[100], line[100], obs_file[100], dump_file[100], com_dump[100];
	erase_all_model();
	FILE *digit_dump;
	int total_files_trained = 20;

	int tsize = 20;

	for(int d = 0; d<=9; d++){
		erase_model();

		sprintf(dump_file, "results/training/training_digit_%d.txt", d);
		FILE *dig_dump = fopen(dump_file, "w");

		fprintf(common_dump, " * DIGIT %d \n", d);
		fprintf(dig_dump, " * DIGIT %d \n", d);

		for(int u = 1; u <= tsize; u++){

			sprintf(filename, "input/dataset/234101006_E_%d_%d.txt", d, u);

			FILE *f = fopen(filename, "r");

			if(f == NULL){
				printf("Issue in opening file %s", filename);
				exit(1);
			}

			fprintf(dig_dump, "\nOpening file %s =============================\n", filename);
			fprintf(common_dump, "\nOpening file %s =============================\n", filename);

			//setting dcshift and nfactor
			setupGlobal(filename);

			sSize = 0;
			//reading the samples and normalizing them
			while(!feof(f)){
				fgets(line, 100, f);

				//input file may contain header, so we skip it
				if(!isalpha(line[0])){
					int y = atof(line);
					double normalizedX = floor((y-dcShift)*nFactor);
					//if(abs(normalizedX) > 1)
					sample[sSize++] = normalizedX;
				}
			}
			fclose(f);

			//framing
			//generating observation seq
			sprintf(obs_file, "output/obs_seq/HMM_OBS_SEQ_%d_%d.txt", d, u);
			get_obs_sequence(obs_file);
			fprintf(dig_dump, "->obs seq: ");
			fprintf(common_dump, "->obs seq: ");

			for(int i=1; i<=T; i++){
				fprintf(dig_dump, "%4d ", O[i]);
				fprintf(common_dump, "%4d ", O[i]);
			}

			fprintf(dig_dump, "\n");
			fprintf(common_dump, "\n");

			//initializing model
			if(train == 0)
				initialize_lambda(d, 1, "--");
			else
				initialize_lambda(d, 1, "avg");

			int iteration = 1;
			//starts converging model upto CONVERGE_ITERATIONS or till convergence whichever reach early
			pstar = 0, prev_p_star = -1;
			while(pstar > prev_p_star && iteration < 1000){
				
				iteration++;
				prev_p_star = pstar;
				forward_procedure();
				backward_procedure();
				viterbi();

				//printing in log file
				fprintf(dig_dump, "iteration: %d\n", iteration);
				fprintf(dig_dump, "-->pstar : %g\n", pstar);
				fprintf(dig_dump, "-->qstar : ");
				for(int i=1; i<=T; i++){
					fprintf(dig_dump, "%d ", Q[i]);
				}
				fprintf(dig_dump, "\n");

				get_xi();
				get_gamma();
				re_estimation();
			}

			//writing final state sequence
			fprintf(common_dump, "-->qstar : ");
			for(int i=1; i<=T; i++){
				fprintf(common_dump, "%d ", Q[i]);
			}
			fprintf(common_dump, "\n");

			//writing final model in the log file
			fprintf(dig_dump, "-------------------------------Final Model Lambda (Pi, A, B) after iterations %d--------------------------------\n", iteration);
			fprintf(common_dump, "-------------------------------Final Model Lambda (Pi, A, B) after iterations %d--------------------------------\n", iteration);
			save_converged_model(dig_dump);
			save_converged_model(common_dump);

			add_to_avg_model();
			save_final_model(u, d);
		}
		fclose(dig_dump);
		average_of_avg_model(tsize);
		save_avg_model(d);
		erase_avg_model();

		//system("pause");
	}
	train++;
}

//To access the codebook from a file.
void read_codebook(){
	ifstream in("codebook1.txt");
	for (int i = 0; i < CB_SIZE; i++){
		for (int j = 0; j < P; j++){
			in >> codeBook[i][j];
		}
	}
	in.close();
}

// Predictions are generated by loading the model and executing the solution for problem 1
void test_prediction(){
	test_ans = 0;
	max_pobs_model = 0;
	for(int k = 0; k<=9; k++){
		read_average_model(k);
		solution_to_prob1(k);
	}

	printf("Digit Recognized =  %d\n", test_ans);
}

//performs live prediction of the data
void real_time_testing(){
	char obs_file[100], line[100];
	printf("\n----------Live testing----------\n");

	system("Recording_Module.exe 3 input.wav input_file.txt");

	initialize_lambda(0, 0, "--");

	FILE *f = fopen("input_file.txt", "r");
	if(f == NULL){
		printf("Issue in opening file input_file.txt");
		exit(1);
	}

	//setting dcshift and nfactor
	setupGlobal("input_file.txt");

	sSize = 0;
	//reading the samples and normalizing them
	while(!feof(f)){
		fgets(line, 100, f);

		//input file may contain header, so we skip it
		if(!isalpha(line[0])){
			int y = atof(line);
			double normalizedX = floor((y-dcShift)*nFactor);
			//if(abs(normalizedX) > 1)
			sample[sSize++] = normalizedX;
		}
	}
	fclose(f);
	get_obs_sequence("output/live_test/obs_seq.txt");

	test_prediction();
}



//function to test the models
void testing(){
	char filename[100], line[100], test_file[100];
	int correctAns = 0, totalCorrectAns = 0;

	int accuracy_array[10] = {0};

	for(int d = 0; d<=9; d++){
		sprintf(test_file, "results/testing/offline/offline_testing_digit_%d.txt", d);
		FILE *fp = fopen(test_file, "w");
		correctAns = 0;
		fprintf(fp, "--------------------------------------------* Digit %d *--------------------------------------------------------\n", d);

		for(int j = 21; j<=30; j++){
			sprintf(filename, "input/dataset/234101006_E_%d_%d.txt", d,j);
			printf("\n\n--------Reading input from the file: %s------\n", filename);

			FILE *f = fopen(filename, "r");
			if(f == NULL){
				printf("Issue in opening file input_file.txt");
				exit(1);
			}
			fprintf(fp, "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n");
			fprintf(fp, "\n---------> Reading file %s <---------\n", filename);

			//setting dcshift and nfactor
			setupGlobal(filename);

			sSize = 0;
			//reading the samples and normalizing them
			while(!feof(f)){
				fgets(line, 100, f);

				//input file may contain header, so we skip it
				if(!isalpha(line[0])){
					int y = atof(line);
					double normalizedX = floor((y-dcShift)*nFactor);
					sample[sSize++] = normalizedX;
				}
			}
			fclose(f);

			//generating observation sequence
			get_obs_sequence("output/live_test/obs_seq.txt");

			fprintf(fp, "observation seq obtained -- ");
			for(int i=1; i<=T; i++){
				fprintf(fp, "%d\t", O[i]);
			}
			fprintf(fp, "\n");

			test_ans = 0;
			max_pobs_model = 0;
			for(int k = 0; k<=9; k++){
				read_average_model(k);
				solution_to_prob1(k, fp);
				erase_avg_model();
			}

			printf("\nDigit Recognized = %d\n", test_ans);
			printf("Original Digit = %d\n", d);

			fprintf(fp, "Digit Recognized: %d\n", test_ans);
			fprintf(fp, "Original Digit: %d\n", d);
			if(test_ans == d) correctAns++, totalCorrectAns++;
		}
		//printf("Accuracy for the digit %d is : %lf % \n", d, (correctAns / 10.0f)*100);
		accuracy_array[d] = correctAns ;
		fprintf(fp, "Accuracy for the digit %d is : %lf % \n", d, (correctAns / 10.0f)*100);
	
		fclose(fp);
	}
	printf("\n\n");
	for(int i = 0;i<=9; i++){
		cout<<"Accuracy for the digit "<<i<<" is "<<(accuracy_array[i]/10.0)*100<<"%"<<endl;
		//printf("Accuracy for the digit %d is : %f % \n", i, (accuracy_array[i] / 10.0)*100);
	}

	cout<<"System Accuracy is "<<totalCorrectAns<<"%\n"<<endl;
	//printf("System Accuracy: %d %\n\n", totalCorrectAns);
}



void menu() {

	printf("\n\n");
    char choice;
    while (true) {
        printf("Enter 1. for automated test on test files\n"
               "Enter 2. for live testing\n"
               "Enter your choice: ");
        scanf(" %c", &choice);

        switch (choice) {
            case '1':
                
                printf("\nAutomated testing...\n");
				testing();
                break;

           
            case '2':
                if (environment_known == 0) {
                    printf("--------------Recording silence--------------\n");
                    system("Recording_Module.exe 3 silence.wav silence_file.txt");
                    environment_known = 1;
                }
                is_live_testing = 1;
                printf("Real-time testing...\n");
				real_time_testing();
                is_live_testing = 0;
                break;

            case '0':
                printf("Exiting the program\n");
                return;

            default:
                printf("Invalid choice. Try again.\n");
        }
    }
}


int _tmain(int argc, _TCHAR* argv[])

{
	char com_dump[100];
	sprintf(com_dump, "results/training/common_dump.txt");
	common_dump = fopen(com_dump, "w");

	read_codebook();
	training();
	menu();
	
	return 0;
}
