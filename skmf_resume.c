/*
Recreation and execution of the proper code's IDEAS coded in my taste. Hopefuly works!
I liberally took some parts that were better if left alone.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>


#define Z 12
#define const_R 8.31446261815324
#define M_RAN_INVM32 2.32830643653869628906e-010


typedef struct sites{
	int x;
	int y;
	int z;
	int r2;
	int pos;
} sites;

typedef struct properties{
	double c;
	double dc;
	double pn_sumc;
	double M;
	double V;
	int numn_temp;
} properties;

int comparator(const void *p1, const void *p2){
	const sites *a=p1;
	const sites *b=p2;
	return ((a->r2)-(b->r2));
};


//STOLEN CONTENT
// data type to store the state of the pseudo RNG
typedef struct xorshift128_state
{
  uint32_t x, y, z, w;
} prng_state;

// pseudo random number generator
/* Algorithm is an implementation of "xor128" from Marsaglia, "Xorshift RNGs" */
/* At initialization the state array should have at least 1 non-zero element */
double xorshift128(struct xorshift128_state *state)
{
	uint32_t t = state->x;
	uint32_t const s = state->w;

	state->x = state->y;
	state->y = state->z;
	state->z = state->w;

	t ^= t << 11;
	t ^= t >> 8;
	state->w = t ^ (s ^ (s >> 19));

	/* Conversion to float is based on:
	Doornik, J.A. Conversion of high-period random numbers to floating point.
	ACM Trans. Model. Comput. Simul. 17, 1, Article 3 */
	return( ((int) state->w) * M_RAN_INVM32 + 0.5);
}


int main(int argc, char* argv[]){

	int save=-1;
	int particles=0;
	int savestate=0;
	int cursorpos=0;
	int cposlag=0;
	char buffer1[50];
	int buffer2;
	char pbuff[50];
	int curstart=0;
	int curend=0;
	int partstart=0;


//variables for savescumming

	double c_start;
	double T;
	double V0, V1, M0, M1;
	double V0a, V0b, V1a, V1b;
	double Eps_0_0, Eps_0_1;
	double Eps_E_0, Eps_E_1;
	double Eps_S_0, Eps_S_1;
	double Sig_E_0, Sig_E_1;
	double Sig_S_0, Sig_S_1;

//variables that will be read.

	double R;
	int Ri, Rs;


	int atomistic_deposition = 1;
	int deposition_on;
	int extra_saves;
	int Ss, Ns;
	int Ds=0;
	double deporate;
	int S=0;
	int Si=0;
	double sumni;
	double diff;
	int str_buf = 500;
	char line[500];
	char name[500];
	double tau=0;
	double dtau;
	double An_tilde;
	int substrate_plus;
	int N=0;
	double Eij;
	int seed;
	prng_state RNGstate;
	int *posindex;
	int *numn;
	int **nn;

	numn = malloc((N+1) * sizeof(int));
	nn = (int**) malloc ((N+1)*sizeof(int*));
	nn[N] = (int*) malloc (Z*sizeof(int));

	sites *sitepos;
	properties *sitecomp;

	sitepos = (sites*) malloc((N+1) * sizeof(sites));
	sitecomp = (properties*) malloc((N+1) * sizeof(properties));

	FILE *input = fopen("input.dat","r");

	if(input == NULL){
		printf("No input file present.\n");
		return -1;
	};

	// reading description line; reading variable name and value
	fgets(line,str_buf,input); fscanf(input,"%s %d\n",name,&seed);

	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&T);

	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&V0a);
	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&V0b);
	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&V1a);
	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&V1b);

	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&M0);
	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&M1);

	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&Eps_0_0);
	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&Eps_0_1);

	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&Eps_E_1);
	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&Eps_S_0);
	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&Eps_S_1);
	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&Eps_E_0);

	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&Sig_E_0);
	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&Sig_E_1);
	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&Sig_S_0);
	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&Sig_S_1);

	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&An_tilde);
	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&dtau);
	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&c_start);

	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&R);
	fgets(line,str_buf,input); fscanf(input,"%s %d\n",name,&Ss);
	fgets(line,str_buf,input); fscanf(input,"%s %d\n",name,&Ns);
	fgets(line,str_buf,input); fscanf(input,"%s %d\n",name,&deposition_on);
	fgets(line,str_buf,input); fscanf(input,"%s %d\n",name,&atomistic_deposition);
	fgets(line,str_buf,input); fscanf(input,"%s %lf\n",name,&deporate);
	fgets(line,str_buf,input); fscanf(input,"%s %d\n",name,&extra_saves);
	fgets(line,str_buf,input); fscanf(input,"%s %d\n",name,&substrate_plus);

	fclose(input);	// closing input file

	//THIS IS WHERE THE FIRST PART OF THE SAVE LOAD FUNCTION RESIDES
	
	if(argc>1){ //if we give input

	FILE *old=fopen("output.xyz","r+");

	if(old==NULL){
		printf("ERROR, no output.xyz file present\n");
		return -2;
	}
		
	
	if(argc>1){
		sscanf(argv[1],"%d",&save);
		printf("Resume at step: %d\n",save); 
	};
	
	if((save!=-1)&&(old!=NULL)){
		
		int i=0;

		while(!feof(old)){ 
				
			for(int j=0;j<50;j++){ //we lag behind to remember the string from last line (matters when the last line is the number of atoms)
				pbuff[j]=buffer1[j];
			}
				
			fscanf(old,"%s %d %*s %*s\n",buffer1,&buffer2); //scan for first two entities in a line
			cposlag=cursorpos; //we lag the cursor position, to remember the beginning of this line too
			cursorpos=ftell(old);  //update the saved cursor position to the start of the next line
			
			
			if(!strcmp("save:",buffer1)){
				cposlag=cposlag-strlen(pbuff)-1; //when we reach the save flag, we move our cursor back the length of the number of atoms + the enter, so we have the cursor position of the start of the nth savefile before number of atoms
				savestate=buffer2;
				particles=strtol(pbuff,NULL,10); //we obtain number of atoms
			}
		
			if(savestate==(save-1)){curstart=cposlag;} //save where we need to start reading next time
			if(savestate==save){curend=cposlag;partstart=particles;} //save where we will stop reading and resume writing, and the number of atoms in the particle
		
			/*if((i<5)){
				printf("buffers:\n%s\n%d\n",buffer1,buffer2);
				printf("particles:%d, savestate:%d, cursorpos:%d\n",particles,savestate,cposlag);
			}//*/	
				
			i++;
		}
	}else{
		printf("ERROR, invalid savestate input\n"); //just filters out not number inputs
		return -3;
	}
		
	fclose(old);
}//if argc	

	printf("start cursor at:%d\n",curstart); //if this and a later cursor position don't match the user knows there is a problem





	printf("substrate plus:%d\n",substrate_plus);
	// I wanted to make the programs compatible, so I copied the names, I rewrote


	V0 = (V0a + V0b*T)/(Z*const_R*T);
	V1 = (V1a + V1b*T)/(Z*const_R*T);
	An_tilde = An_tilde * sqrt(3.) / sqrt(dtau);
double depodelay= 1.0/dtau/deporate;

printf("depodelay:%g\n",depodelay );

	if(deporate > (1.0/dtau)){
		printf("Insufficient deposited atoms for each timestep\n");
		return -1;
	};

	Ri=(int) ceil(R);

	printf("Ri:%d\n",Ri);

	//boxing the particle
		int boxwidth = 2*Ri+1;
		int boxheight = Ri+1;
		int blayer = boxwidth*boxwidth;
		int bvolume = boxwidth*boxwidth*boxheight;

		posindex = malloc((bvolume)*sizeof(int));

		for(int i=0;i<bvolume;i++){
			posindex[i]= INT_MAX;
		};

	int r_neigh_indices[Z];
	int counter=0;
	for(int i=-1;i<2;i++){
		for(int j=-1;j<2;j++){
			for(int k=-1;k<2;k++){
			if((abs(i)+abs(j)+abs(k))==2){
			r_neigh_indices[counter]= i + j*blayer + k*boxwidth;
			counter++;
			};//if
		}; //for k
	}; // for j
}; // for i
counter=0;
//nem ugyan az a szomszéd index sorrendje mint ami az eredeti kódban



//az rng-nek utána olvastam, tök érdekes volt, ki is próbáltam egy megírt játékban kicseréltem a rand függvényt erre, ott működött, ezért innen meg átemeltem mert nem kell nekem belenyúlkálni ha működik. Meg amúgy is standard implementációs módja van, ezt ha átírom nem lesz jobb semmi.

if (seed == 0){
	srand(time(NULL));
	RNGstate = (prng_state){.x=(uint32_t) rand(),.y=(uint32_t) rand(),.z=(uint32_t) rand(),.w=(uint32_t) rand()};
}else{
	RNGstate = (prng_state){.x=12345,.y=246810,.z=3691215,.w=48121620};
};


//Classic way to initialize used in every code thus far.

for(int i=-Ri;i<=Ri;i++){
	for(int j=-Ri;j<=Ri;j++){
		for(int k=0;k<=Ri;k++){
			int distance = i*i+k*k+j*j;

			if(distance< (R*R)){
				if((abs(i)+abs(j)+abs(k))%2==0){

					sitepos = (sites*) realloc (sitepos, (N+1)*sizeof(sites));
					numn = (int*) realloc (numn, (N+1)*sizeof(int));
					nn= (int**) realloc (nn, (N+1)*sizeof(*nn));
					nn[N]= (int*) malloc (Z*sizeof(int));

					sitepos[N].x=i;
					sitepos[N].y=j;
					sitepos[N].z=k;
					sitepos[N].r2=distance;
					sitepos[N].pos= i+j*boxwidth+k*blayer+Ri*boxwidth+Ri;

					numn[N]=0;

					for(int l=0;l<Z;l++){
					nn[N][l]=INT_MAX;

				  }; //neighbours

					N++;
				};//is the site on the FCC lattice?
			};//is the site in the dome's range? (is it in the particle?)

		};//for k
	};//for j
};//for i

printf("%d atoms are present\n",N);


qsort(sitepos,N,sizeof(sites),comparator);


sitecomp= (properties*) realloc (sitecomp,N*sizeof(properties));


for(int i=0;i<N;i++){
	posindex[sitepos[i].pos]=i;

	if(atomistic_deposition==1){
		int atoms= (int) floor(c_start*N);
		if(i<atoms)					{sitecomp[i].c=1.0;}
		else if(i>atoms)		{sitecomp[i].c=0.0;}
		else								{sitecomp[i].c=c_start*N -atoms;};
	}
	else{
		sitecomp[i].c=c_start;
	};
	sitecomp[i].dc = 0.0;
	sitecomp[i].pn_sumc = 0.0;
	sitecomp[i].M = 0.0;
	sitecomp[i].V = 0.0;
	sitecomp[i].numn_temp = 0;
};//for i


printf("deposition initialized\n");

if(atomistic_deposition==1){
	for (int i=N-1;i>0;i--){
		int j = (int) floor(xorshift128(&RNGstate) * i);
		double temp = sitecomp[i].c;
		sitecomp[i].c = sitecomp[j].c;
		sitecomp[j].c = temp;
	}; //for i
}; //if


for (int i=0; i<N; i++){
	for (int k=0; k<Z; k++){
		int npos = sitepos[i].pos + r_neigh_indices[k];
		int j;
		if (0 <= npos && npos < bvolume){
			j = posindex[npos];
		}else {j = INT_MAX;}; // end of else
		if (j < INT_MAX){
			nn[j][numn[j]] = i; // adding site i to j's neighbours list
			numn[j]++; // increment number of neighbours of site j
		};// if j<int_max
	}; //for k

	printf("\rNeighbouring:%.1f%%",100.*i/(N-1));
	fflush(stdout);
}; // for i



printf("\n");
free(posindex);
int depositioncounter=0;
if (deposition_on == 1){
	printf("%d time step per atom\n",(int) ceil(depodelay));
}else{
	depositioncounter=N;
	for (int i=0;i<N;i++){
		sitecomp[i].numn_temp = numn[i];
	};//for
	printf("initially complete particle\n");
};//else


//SAVE RESUME PART 2: ELECTRIC BOOGALOO

char* flag[5];

if(argc>1){*flag="r+";}else{*flag="w";} //if there is an input to the file we continue, if there isn't we overwrite the old xyz file



FILE *fout=fopen("output.xyz",*flag); 


	if(argc>1){ //only do this wrap around bs if necessary
	
		double conc;
		int counttodepo=0;
		int skip=0;
	
		fseek(fout,curstart,SEEK_SET);
	
		printf("Successfuly found start position: %ld\n",ftell(fout)); 
		printf("Loading %d atoms\n",partstart);
	
		S=save;
		depositioncounter=partstart-(((2*Ri+2*substrate_plus+1)*(2*Ri+2*substrate_plus+1))/2); //I don't deposit substrate
		
		while(ftell(fout)<curend-1){
	
			if(skip>1){
				fscanf(fout,"%*d %*d %*d %lf\n",&conc);
				 
				if(conc!=-1){
				
				sitecomp[counttodepo].c=conc;   //we read the concentration from the file
				
				//printf("numn_count:%d\n",numn[counttodepo]);
				
				//we have to establish the neighbour relation of already deposited atoms
				for(int i=0;i<numn[counttodepo];i++){
				sitecomp[nn[counttodepo][i]].numn_temp++;
				}//for i */
				
				counttodepo++;
				}
			
			
			}else{fgets(line,str_buf,fout);}	
			skip++;
			if(counttodepo==depositioncounter){break;} //to make sure we break while, for some reason the ftell statement worked in my separate sandbox code, but not in the live skmf one 
		}//while*/
	
		printf("cursor ended up at:%ld\n",ftell(fout));
		printf("cursor should be at:%d\n",curend);
	
	}
	
	// open output file

    

	//i don't know how to write ovito files, I read the documentation and used the already existing programs as a base CORRECTION NOW I KNOW HOW TO WRITE OVITO FILES


	while (S <= Ns){
		if ((Si==Ss)) {

			fprintf(fout,"%d\n",(depositioncounter+(((2*Ri+2*substrate_plus+1)*(2*Ri+2*substrate_plus+1))/2)));	// number of particles, substrate magic included
			
			fprintf(fout,"save: %d ",S);
			
			fprintf(fout,"Lattice=\"");
			fprintf(fout,"%d %d %d " ,0,0,0);		// x_vector of bounding box
			fprintf(fout,"%d %d %d " ,0,0,0);	// y_vector of bounding box
			fprintf(fout,"%d %d %d\"",0,0,0);		// z_vector of bounding box

			fprintf(fout,"Properties=pos:I:3:Color:R:1 ");
			fprintf(fout,"Time=%lf\n",tau);


			for(int i=-(Ri+substrate_plus); i <= (Ri+substrate_plus); i++){
				for(int j=-(Ri+substrate_plus); j <= (Ri+substrate_plus); j++){
					if((i+j-1)%2==0){
						fprintf(fout, "%d %d -1 -1\n",i,j);
					};//if
				};//for j
			};//for i

			for(int i=0;i<depositioncounter;i++){
					fprintf(fout,"%d %d %d %lf\n", sitepos[i].x,sitepos[i].y,sitepos[i].z,sitecomp[i].c);
				};//for i

				if(deposition_on == 1){
					if (depositioncounter < N){
						printf("\rDeposition underway %d of %d records",S,Ns);
						fflush(stdout);
					}else{
						printf("\rFinal particle simulation %d of %d records",S,Ns);
						fflush(stdout);
					};//inner if's fi
				}else{
					printf("\rRunning %d of %d records",S,Ns); //NOT FINAL
					fflush(stdout);
				};//large if

			S++;		// increase save counter
			Si = 0;	// set cycle counter to 0
		};
		// actual saving


		if(deposition_on==1){
			if((depodelay<=Ds) && (depositioncounter<N) || (depositioncounter==0)){
				for(int i=0;i<numn[depositioncounter];i++){
					sitecomp[nn[depositioncounter][i]].numn_temp++;
                    Si=1;
				};// for i
				Ds=0;
				depositioncounter++;

			}//if depodepo
			else if( (depodelay<=Ds)  && (depositioncounter>=N)){
				extra_saves--;
				if(extra_saves<=0){
					printf("simulation over\n");
					return 0;
				};//if extra saves
				Si=1;
				Ds=0;
			};//else
		} 
		else{
			Si=1;			
		};//fi deposition on;
		//neighbouring

		for (int i=0; i<depositioncounter; i++){
			sitecomp[i].pn_sumc = 0.0;
			for (int j=0; j<sitecomp[i].numn_temp; j++){
				int in = nn[i][j];
				//if(in!=INT_MAX){
				sitecomp[i].pn_sumc += sitecomp[in].c;
				//};
			};//for j
			sitecomp[i].V = V0 + V1 * (2.*sitecomp[i].c - 1.);
			sitecomp[i].M = M0 + M1 * (2.*sitecomp[i].c - 1.);
		};//for i


		//stroking the new equations

		for(int i=0; i<depositioncounter; i++){
			double c_i=sitecomp[i].c;
			int np_i=sitecomp[i].numn_temp;
			int ns_i=0;
			if(sitepos[i].z==0){ns_i=4;};
			int ne_i=Z-np_i-ns_i;

			for(int in=0; in<sitecomp[i].numn_temp; in++){
				int j=nn[i][in];

				double c_j=sitecomp[j].c;
				double np_j=sitecomp[j].numn_temp;
				int ns_j=0;
				if(sitepos[j].z==0){ns_j=4;};
				int ne_j=Z-ns_j-np_j;

				Eij = (sitecomp[i].M - sitecomp[i].V)*sitecomp[i].pn_sumc + (sitecomp[j].M + sitecomp[j].V)*sitecomp[j].pn_sumc;

				Eij += Z*(2.*Eps_0_0 + 2.*Eps_0_1*(c_i + c_j - 1.) + M1*(c_i - c_j));

				Eij += (1./2.)*V0*(np_i - np_j) + (1./2.)*V1*(np_i*(2.*c_i - 1.) - np_j*(2.*c_j - 1.));

				Eij += Eps_E_0*(ne_i + ne_j) + Eps_E_1*(ne_i*(2.*c_i - 1.) + ne_j*(2.*c_j - 1.));
				Eij += Eps_S_0*(ns_i + ns_j) + Eps_S_1*(ns_i*(2.*c_i - 1.) + ns_j*(2.*c_j - 1.));
				Eij += Sig_E_0*(ne_i - ne_j) + Sig_E_1*(ne_i*(2.*c_i - 1.) - ne_j*(2.*c_j - 1.));
				Eij += Sig_S_0*(ns_i - ns_j) + Sig_S_1*(ns_i*(2.*c_i - 1.) - ns_j*(2.*c_j - 1.));

				double dci = - c_i * (1.0 - c_j) * exp(Eij) * dtau;
				sitecomp[i].dc += dci;
				sitecomp[j].dc -= dci;
		//		};//if j has been deposited
			};//for in
		};//big for i

		//is the cycle stable?

		for (int i=0; i<depositioncounter; i++){
			if ( ((sitecomp[i].c + sitecomp[i].dc) < 0.0) || (1.0 < (sitecomp[i].c + sitecomp[i].dc)) ){
				printf("\nKMF cycle unstable\n");
				printf("Decrease dtau.\n");
				return -1;
			};
		};

		//we put the S in SKMF

		if (An_tilde > 0.){
			for (int i=0; i<depositioncounter; i++){
				for (int in=0; in<sitecomp[i].numn_temp; in++){
					double rn = xorshift128(&RNGstate);
					int j = nn[i][in];
					double dci = - sitecomp[i].c * (1. - sitecomp[j].c) * An_tilde * (2.*rn - 1.) * dtau;
					sitecomp[i].dc += dci;
					sitecomp[j].dc -= dci;
				};
			};
		};

		//do the composition shifting

		for(int i=0; i<depositioncounter; i++){
			sitecomp[i].c += sitecomp[i].dc;
			sitecomp[i].dc=0.0;
		};


		//redistribution

		if (An_tilde > 0.0){
			for (int i=0; i<depositioncounter; i++){
				double diff;
				if (sitecomp[i].c < 0.0){
					printf("\nlow\n");
					double sumni = 0.0;

					for (int j=0;j<sitecomp[i].numn_temp;j++){
						sumni += sitecomp[nn[i][j]].c;
					};//for j

					diff = sitecomp[i].c;
					sitecomp[i].c = 0.;

					for (int j=0; j<sitecomp[i].numn_temp; j++){
						double dci= diff * sitecomp[nn[i][j]].c/sumni;
						sitecomp[nn[i][j]].c+= dci;
					}// for j
				}// if sitecomp[i]
				else if (1. < sitecomp[i].c){
				printf("\nhigh\n");
					double sumni = 0.0;
					for (int j=0; j<numn[i]; j++){
						sumni+= sitecomp[nn[i][j]].c;
					}//for j
					diff= sitecomp[i].c-1.0;
					sitecomp[i].c = 1.0;

					for (int j=0; j<sitecomp[i].numn_temp; j++){
						double dci= diff * (1. - sitecomp[nn[i][j]].c ) / (sitecomp[i].numn_temp - sumni);
						sitecomp[nn[i][j]].c += dci;
					}// for j
				}//else if
			}//for i


			for(int i=0; i<depositioncounter; i++){
				if ( (sitecomp[i].c < 0. ) || ( 1. < sitecomp[i].c )){
					printf("Composition is out of range after backdistributing.\n");
					printf("Decrease the value of dtau.\n");
					return -2;
				}
				if (sitecomp[i].c != sitecomp[i].c ){
					printf("Composition is not a number.\n");
					printf("Probably one of the input parameters is not a number.\n");
				}
			}//for i


		};// if an_tilde

		Ds++;

		tau+=dtau;
};// save while

fclose(fout);
printf("\nSimulation complete\n");

return 0;
};
