//Calculating the DC shift.
void get_DC_shift(){
	long int sample_count = 0;
	int cnt = 0;
    FILE *fp;
    char line[80];
	double cValue;
	double cEnergy = 0;

    //reading dc_shift.txt file
	if(is_live_testing == 0)
		fp = fopen("silence.txt", "r");
	else
		fp = fopen("silence_file.txt", "r");

    if(fp == NULL){
        printf("Silence File not found\n");
        exit(1);
    }

	dcShift = 0;
	silenceEnergy = 0; //resetting the silence Energy to 0
    while(!feof(fp)){
        fgets(line, 80, fp);
		cValue = atof(line);
        dcShift += cValue;

		if(cnt == 100){
			if(silenceEnergy < cEnergy){
				silenceEnergy = cEnergy;
			}
			cEnergy = 0;
			cnt = 0;
		}
		cnt++;
		cEnergy += cValue * cValue;

        sample_count++;
    }
    dcShift /= sample_count;

    fclose(fp);

}

//A function to initialize global variables such as `max` and `nFactor` based on the characteristics of the vowel recording file.
//`max` and `nFactor` are determined based on the vowel recording file and are employed for normalization purposes.
void setupGlobal(char *filename){
    FILE *fp;
    long int totalSample = 0;
    char line[100];

    fp = fopen(filename, "r");
    if(fp == NULL){
        printf("File not found\n");
    }

    //Retrieve the maximum value.
    mx = 0;
    while(!feof(fp)){
        fgets(line, 100, fp);
        if(!isalpha(line[0])){
            totalSample++;
            if(mx < abs(atoi(line)))
                mx = abs(atoi(line));
        }
    }

    nFactor = (double)LIMIT/mx;

    //Establishing the DC shift.
    get_DC_shift();
    fclose(fp);
}

//Computing Tokhura's distance using the codebook.
void calculate_tokhura_distance(long double cepstralCoeff[12], int index, FILE *fp) {
    int min_index = 0;
    long double min_distance = DBL_MAX;
    
    for (int j = 0; j < CB_SIZE; j++) {
        long double distance = 0;
        for (int i = 0; i < P; i++) {
            distance += tokhuraWeights[i] * pow(cepstralCoeff[i] - codeBook[j][i], 2);
        }

        if (distance < min_distance) {
            min_distance = distance;
            min_index = j;
        }
    }

    O[index] = min_index + 1;

    fprintf(fp, "%4d ", O[index]);
}


//This function calculates the cepstral coefficients (Ci's).
void get_Ci(){
	//if(fake == 62) system("pause");
	double sum=0;
	Ci[0]=log(Ri[0]*Ri[0]);

	for(int m=1;m<=P;m++){
		sum=0;
		for(int k=1;k<m;k++){
			sum += (k*Ci[k]*Ai[m-k])/(m*1.0);
		}
		Ci[m]=Ai[m]+sum;
	}
	//fake++;
}

// A function to apply the Durbin algorithm and determine the values of the coefficients (ai's).
void run_durbin_algo(){
	double alpha[13][13],E[13],K[13];
	double sum = 0;
	E[0] = Ri[0];
	//loop for p from 1 to 12
	for(int i=1;i<=P;i++){
		sum=0;
		for(int j=1;j<=i-1;j++){
			sum += alpha[i-1][j]*Ri[i-j];
		}

		K[i]=(Ri[i]-sum)/E[i-1];

		alpha[i][i]=K[i];

		for(int j=1;j<=i-1;j++){
			alpha[i][j]=alpha[i-1][j] - K[i]*alpha[i-1][i-j];
		}

		E[i]=(1-(K[i]*K[i]))*E[i-1];
	}

	//storing the ai values
	for(int i=1;i<=P;i++){
		Ai[i]= alpha[P][i];
	}
}

//Computing the values of Ri.
void get_Ri(double *samp){
	//if(fake == 62) system("pause");
	for(int m =0; m<=P; m++){
		Ri[m] = 0;
		for(int k=0; k<FRAME_SIZE-m; k++){
			Ri[m] += samp[k]*samp[k+m];
		}
	}
}

//A function to implement the Raised Sine Window on \(C_i\) for each frame.
void RaisedSineWindow(){
	long double sum=0;
	for(int m=1;m<=P;m++){
		sum = (P/2)*sin((PI*m)/P);
		Ci[m]*=sum;
	}
}

//Computing the values for c'
void get_C_dash(double *samp){
	get_Ri(samp);
	run_durbin_algo();
	get_Ci();
	RaisedSineWindow();

}

void trim_digit_wave(){
	int num_frames = 0;
	int cnt =0;
	enSize = 0;
	double cEnergySum = 0, multiplier = 3;
	int startMarker = -1, endMarker = -1;

	for(int i=0; i<sSize; i++){
		double cEnergy = sample[i]*sample[i];
		if(cnt == 100){
			energy[enSize++] = cEnergySum;
			cEnergySum = 0;
			cnt = 0;
		}
		cnt++;
		cEnergySum += cEnergy;
	}

	
	int min_samples = 11200;

	for(int i=0; i<enSize-4; i++){
		if(startMarker == -1 && endMarker == -1 && energy[i+1] > multiplier * silenceEnergy && energy[i+2] > multiplier * silenceEnergy && energy[i+3] > multiplier * silenceEnergy && energy[i+4] > multiplier * silenceEnergy){
			startMarker = i*100;
		}else if(startMarker != -1 && endMarker == -1 && energy[i+1] <= multiplier * silenceEnergy && energy[i+2] <= multiplier * silenceEnergy && energy[i+3] <= multiplier * silenceEnergy && energy[i+4] <= multiplier * silenceEnergy){
			int diff = i*100 - startMarker;
			if(diff < min_samples){
				startMarker = 0 > (startMarker - (min_samples - diff)/2) ? 0 : (startMarker - (min_samples - diff)/2);
				endMarker = enSize*100 < (i*100 + (min_samples - diff)/2) ? enSize*100 : (i*100 + (min_samples - diff)/2);
			}
			else
				endMarker = i*100;
		}else if(startMarker != -1 && endMarker!= -1) break;
	}

	sampSize = 0;
	ofstream out("trim.txt");
	for(int i=startMarker; i<=endMarker; i++){
		steadySamp[sampSize++] = sample[i];
		out<<sample[i]<<endl;
	}
	out.close();
	//system("pause");
}


//Create a sequence of observations.
void get_obs_sequence(char *filename){
	int obs_ind = 1;
	FILE *op = fopen(filename, "w");
	if(op == NULL) {
		printf("File not found\n");
		exit(1);
	}

	trim_digit_wave();
	double fsamp[FRAME_SIZE];
	int num_frames = 0;
	for(int i=0; i<sampSize; i++){
		num_frames++;
		for(int j = 0; j<320; j++)
			fsamp[j] = steadySamp[i++];

		get_C_dash(fsamp);
		calculate_tokhura_distance(Ci, obs_ind++, op);
	}
	T = num_frames;
	
	fprintf(op, "\n");
	fclose(op);
	//cout<<"wrote observation seq in file: "<<filename<<"\n";
}

