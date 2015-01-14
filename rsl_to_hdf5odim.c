/*
	|||||   ||||||||||||||||| SIMEPAR
	|||||   ||||||||||||||||| Sistema Meteorológico do Paraná
	|||||               |||||
	|||||||||||||||||   ||||| Technology and Environmental Information


	(C) Copyleft Sistema Meteorológico Simepar (SIMEPAR)
        http://www.simepar.br
 
	Author: Anderson Luis Gama
	Manager: Cesar Beneti

	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 3 of the License (29 June 
        2007), or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program (LICENCE.txt); if not, write to the Free 
	Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
/*
OPERA Data Information Model (ODIM) format specifications according to http://www.eumetnet.eu/sites/default/files/OPERA2014_O4_ODIM_H5-v2.2.pdf 
*/


#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#define USE_RSL_VARS
#include <rsl.h>
#include <hdf5.h>
#include "NetCDF_CfRadial_int08.h"

//#include "../include/radarreader.h"




/**************** PROTOTYPE ****************/
void radar_to_hdf5odim(Radar* radar, char *outfile);// This is the only function intended to be call externally, it receives a RSL/TRMM radar structure (trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl/) and write a HDF5 file in the ODIM Convention named outfile. In respect to the convention this function expect all volumes of the radar strucure to have the same number of sweeps. Furthermore although compatible with RSL V1.44, only volumes with index smaller than 20 are converted, this is dow to the fact that we are not sure how to name the other volumes in the HDF5 file. To extend this functionaly update the variable rsl_to_name_hdf5

//This functions automate the tedious process of writing to a HDF5 file
void write_attr_text(hid_t loc_id,char *name,char * value);//write atribute text to HDF5
void write_attr_float(hid_t loc_id,char *name,float value);//write atribute float to HDF5
void write_attr_uint(hid_t loc_id,char *name,int value);//write atribute unsignet int to HDF5
void write_attr_double(hid_t loc_id,char *name,double value);//Write atribute double to HDF5
void write_attr_long(hid_t loc_id,char *name,long int value);//write atribute long int to HDF5 64-bits


unsigned char *make_buffer(int moment,Sweep * sweep,double max, double min);//alloc and fill rays*bins buffer with moments data in scaled char to write to HDF5 /datasetXX/moment_YY dataset,  make sure to free it after

int get_unfolding(Ray * ray);//Use prf and prf2 from ray header to calculate unfolding code according with HDF5 convertion
unsigned long long int unix_time(Ray_header header);//convert time in a ray to unix_time (in microseconds)
double get_max(Sweep *sweep);// get the biggers value in sweep;
double get_min(Sweep *sweep);// get the smallest value in sweep;
void get_index_first_ray_in_sweep(Sweep *sweep,int *ray_index); // get the index of the ray with smallest time in sweep;
void get_index_last_ray_in_sweep(Sweep *sweep,int *ray_index); // get the index of the ray with smallest time in sweep;
void get_index_first_ray_in_volume(Volume *volume,int *sweep_index,int *ray_index);// get the index of the ray with smallest time in volume;
void get_index_last_ray_in_volume(Volume *volume,int *sweep_index,int *ray_index);// get the index of the ray with smallest time in volume;

/**************** CODE *********************/

void radar_to_hdf5odim(Radar* radar, char *outfile){
	extern int radar_verbose_flag; /* Defined by RSL */  
	//first of all I need to sort

	RSL_sort_radar (radar);

	char string[100]; //all uses string
	char moment[20],dataset[20];//moment and dataset name string
	hid_t file_id,what_id,how_id,where_id,dataset_id;//specific id's, although dataset_id may change (roll) during the program 
	hid_t group_id, attr_id, dataspace_id,dataarray_id,data_id; //all uses id
	hsize_t dims[2];//dimensions
	int i,j,k,len;// len: length of strings
	int nsweeps, nvolume=0;
	int moments[200];//moment that will be writen, in rsl numeration
	#define MAX_NVOLUME 20
	char *rsl_to_name_hdf5[MAX_NVOLUME]={"DBZH","VRADH","WRADH","DBZH","TH","ZDR","LDR","ZDR","SIGPOW","RHOHV","PHIDP","XZ","ZDR","DBZH","ZDR","DBZH","VRADDH","KDP","TIME","QIND"};
	herr_t status;
	unsigned char * buffer;
	double min,max;
	int temp_int,ray_idx,sweep_idx;
	j=0;
	for(i=0;i<MAX_NVOLUME;i++)if(radar->v[i]!=NULL){moments[j]=i;j++;}//get list of not NULL volumes in Radar
	nvolume=j;//number of volumes to be writen
	nsweeps=radar->v[moments[0]]->h.nsweeps;//number of datasets


	//Overwrite file if already exist
	file_id = H5Fcreate(outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if(file_id<0){ printf("Falling Creating a HDF5 file\n"); exit(-1);}

	//root group
	write_attr_text(file_id,"Conventions","ODIM_H5/V2_2");
	//creat groups
	what_id=H5Gcreate2(file_id,"/what", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//create, white and close atributes of /what
		if(radar->v[moments[0]]->h.nsweeps==1)
			write_attr_text(what_id,"object","SCAN");
		else
			write_attr_text(what_id,"object","PVOL");
		write_attr_text(what_id,"version","H5rad 2.2");
		get_index_first_ray_in_volume(radar->v[moments[0]],&sweep_idx,&ray_idx);
		sprintf (string, "%04i%02i%02i",radar->v[moments[0]]->sweep[sweep_idx]->ray[ray_idx]->h.year,radar->v[moments[0]]->sweep[sweep_idx]->ray[ray_idx]->h.month,radar->v[moments[0]]->sweep[sweep_idx]->ray[ray_idx]->h.day);
		write_attr_text(what_id,"date",string);
		temp_int=radar->v[moments[0]]->sweep[sweep_idx]->ray[ray_idx]->h.sec;
		sprintf (string, "%02i%02i%02i",radar->v[moments[0]]->sweep[sweep_idx]->ray[ray_idx]->h.hour,radar->v[moments[0]]->sweep[sweep_idx]->ray[ray_idx]->h.minute,temp_int);
		write_attr_text(what_id,"time",string);
		sprintf (string, "CMT:%s",radar->h.name);
		write_attr_text(what_id,"source",string);
		

		H5Gclose(what_id);



	where_id=H5Gcreate2(file_id,"/where", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//create, white and close atributes of /where
		write_attr_double(where_id,"lat",radar->h.latd+(radar->h.latm/60.)+(radar->h.lats/3600.));
		write_attr_double(where_id,"lon",radar->h.lond+(radar->h.lonm/60.)+(radar->h.lons/3600.));
		write_attr_double(where_id,"height",radar->h.height);
		H5Gclose(where_id);


	how_id=H5Gcreate2(file_id,"/how", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//create, white and close atributes of /how
		write_attr_double(how_id,"beamwH",radar->v[moments[0]]->sweep[0]->h.vert_half_bw*2);
		write_attr_double(how_id,"beamwV",radar->v[moments[0]]->sweep[0]->h.horz_half_bw*2);
		write_attr_double(how_id,"wavelength",radar->v[moments[0]]->sweep[0]->ray[0]->h.wavelength*100);
		write_attr_long(how_id,"scan_count",nsweeps); 
		H5Gclose(how_id);
		
		

	for(i=0;i<nsweeps;i++){
		int nrays=radar->v[moments[0]]->sweep[i]->h.nrays;
		if(radar_verbose_flag) printf("Writing /dataset%i...\n",i+1);
		sprintf(dataset,"/dataset%d",i+1); 
		dataset_id=H5Gcreate2(file_id,dataset, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		group_id=H5Gcreate2(dataset_id,"what", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		//create, white and close atributes of /datasetXX/what
			write_attr_text(group_id,"product","SCAN");
// not here!!!
			get_index_first_ray_in_sweep(radar->v[moments[0]]->sweep[i],&ray_idx);printf("first ray %i\n",ray_idx);
			sprintf (string, "%04i%02i%02i",radar->v[moments[0]]->sweep[i]->ray[ray_idx]->h.year,radar->v[moments[0]]->sweep[i]->ray[ray_idx]->h.month,radar->v[moments[0]]->sweep[i]->ray[ray_idx]->h.day);
			write_attr_text(group_id,"startdate",string);
			temp_int=radar->v[moments[0]]->sweep[i]->ray[ray_idx]->h.sec;
			sprintf (string, "%02i%02i%02i",radar->v[moments[0]]->sweep[i]->ray[ray_idx]->h.hour,radar->v[moments[0]]->sweep[i]->ray[ray_idx]->h.minute,temp_int);
			write_attr_text(group_id,"starttime",string);
			get_index_last_ray_in_sweep(radar->v[moments[0]]->sweep[i],&ray_idx);printf("last ray %i\n",ray_idx);
			sprintf (string, "%04i%02i%02i",radar->v[moments[0]]->sweep[i]->ray[ray_idx]->h.year,radar->v[moments[0]]->sweep[i]->ray[ray_idx]->h.month,radar->v[moments[0]]->sweep[i]->ray[ray_idx]->h.day);
			write_attr_text(group_id,"enddate",string);
			temp_int=radar->v[moments[0]]->sweep[i]->ray[ray_idx]->h.sec;
			sprintf (string, "%02i%02i%02i",radar->v[moments[0]]->sweep[i]->ray[ray_idx]->h.hour,radar->v[moments[0]]->sweep[i]->ray[ray_idx]->h.minute,temp_int);
			write_attr_text(group_id,"endtime",string);
			//write_attr_uint(group_id,"descriptor_count",nvolume); oh my god where does this go!!!
			H5Gclose(group_id);

		group_id=H5Gcreate2(dataset_id,"where", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		//create, white and close atributes of /datasetXX/where
			write_attr_double(group_id,"elangle",radar->v[moments[0]]->sweep[i]->ray[0]->h.elev);
			write_attr_long(group_id,"nbins",radar->v[moments[0]]->sweep[i]->ray[0]->h.nbins);
			write_attr_double(group_id,"rstart",radar->v[moments[0]]->sweep[i]->ray[0]->h.range_bin1/1000.);
			write_attr_double(group_id,"rscale",radar->v[moments[0]]->sweep[i]->ray[0]->h.gate_size);
			write_attr_long(group_id,"nrays",radar->v[moments[0]]->sweep[i]->h.nrays);
			get_index_first_ray_in_sweep(radar->v[moments[0]]->sweep[i],&ray_idx);
			write_attr_long(group_id,"a1gate",ray_idx);
			write_attr_double(group_id,"startaz",radar->v[moments[0]]->sweep[i]->ray[0]->h.azimuth);
			write_attr_double(group_id,"stopaz",radar->v[moments[0]]->sweep[i]->ray[radar->v[moments[0]]->sweep[i]->h.nrays-1]->h.azimuth);
			printf("startaz %f, stopaz %f\n",radar->v[moments[0]]->sweep[i]->ray[0]->h.azimuth,radar->v[moments[0]]->sweep[i]->ray[radar->v[moments[0]]->sweep[i]->h.nrays-1]->h.azimuth);
			H5Gclose(group_id);

		group_id=H5Gcreate2(dataset_id,"how", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		//create, white and close atributes of /datasetXX/how
			write_attr_long(group_id,"scan_index",i);
			write_attr_long(group_id,"ray_count",radar->v[moments[0]]->sweep[i]->h.nrays);
			write_attr_double(group_id,"lowprf",radar->v[moments[0]]->sweep[i]->ray[0]->h.prf2);
			write_attr_double(group_id,"highprf",radar->v[moments[0]]->sweep[i]->ray[0]->h.prf);
			write_attr_double(group_id,"RXbandwidth",radar->v[moments[0]]->sweep[i]->ray[0]->h.frequency*1000);
			write_attr_double(group_id,"pulsewidth",radar->v[moments[0]]->sweep[i]->ray[0]->h.pulse_width);
			write_attr_double(group_id,"rpm",radar->v[moments[0]]->sweep[i]->ray[0]->h.sweep_rate);
			write_attr_double(group_id,"NI",radar->v[moments[0]]->sweep[i]->ray[0]->h.nyq_vel);
			// Now the azimuth, this is an array, so I must do it by hand
				hid_t attr_dataspace_id,attr_id;
				const hsize_t attr_len=radar->v[moments[0]]->sweep[i]->h.nrays;
				double azimuth[attr_len];

				for(k=0;k<attr_len;k++) azimuth[k]=radar->v[moments[0]]->sweep[i]->ray[k]->h.azimuth-radar->v[moments[0]]->sweep[i]->h.horz_half_bw;
				attr_dataspace_id = H5Screate_simple(1, &attr_len, NULL);
				attr_id=H5Acreate2(group_id,"startazA",H5T_IEEE_F64LE,attr_dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
				H5Awrite(attr_id,H5T_NATIVE_DOUBLE, &azimuth);
				H5Aclose(attr_id);
				H5Sclose(attr_dataspace_id); 

				for(k=0;k<attr_len;k++) azimuth[k]=radar->v[moments[0]]->sweep[i]->ray[k]->h.azimuth+radar->v[moments[0]]->sweep[i]->h.horz_half_bw;
				attr_dataspace_id = H5Screate_simple(1, &attr_len, NULL);
				attr_id=H5Acreate2(group_id,"stopazA",H5T_IEEE_F64LE,attr_dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
				H5Awrite(attr_id,H5T_NATIVE_DOUBLE, &azimuth);
				H5Aclose(attr_id);
				H5Sclose(attr_dataspace_id); 
				

			H5Gclose(group_id);
		for(j=0;j<nvolume;j++){
			if(radar_verbose_flag)printf("        /dataset%i/data%i:  %s\n",i+1,j+1,rsl_to_name_hdf5[moments[j]]);
			sprintf(moment, "data%d", j+1);
			//printf("moment %i",moments[j]);
			//create dataset /datasetXX/moment_YY
			dims[0]=radar->v[moments[j]]->sweep[i]->h.nrays;
			dims[1]=radar->v[moments[j]]->sweep[i]->ray[0]->h.nbins;
			data_id=H5Gcreate2(dataset_id,moment, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			dataspace_id = H5Screate_simple(2, dims, NULL);
			dataarray_id = H5Dcreate2(data_id, "data", H5T_STD_U8LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			//write attribute to /datasetXX/scanYY/what
			max=get_max(radar->v[moments[j]]->sweep[i]);
			min=get_min(radar->v[moments[j]]->sweep[i]);
			group_id=H5Gcreate2(data_id,"what", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			//create, white and close atributes of /datasetXX/dataYY/what
				write_attr_text(group_id,"quantity",rsl_to_name_hdf5[moments[j]]);//moment name
				write_attr_double(group_id,"gain",(max-min)/255.);
				write_attr_double(group_id,"offset",-(max-min)/255.+min);
				write_attr_double(group_id,"nodata",0);
				write_attr_double(group_id,"undetect",0);
			H5Gclose(group_id);
			//get date from sweep, write and close
			buffer=make_buffer(moments[j],radar->v[moments[j]]->sweep[i],max,min);
			status = H5Dwrite(dataarray_id, H5T_STD_U8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
			status = H5Dclose(dataarray_id);
			status = H5Sclose(dataspace_id); 
			H5Gclose(data_id);
			free(buffer);

		}
		status = H5Gclose(dataset_id);
	}
	if(radar_verbose_flag)printf("Closing file...");
	status = H5Fclose(file_id);
	if(radar_verbose_flag)printf("OK\n");
}	


/********* ATTRIBUTE WRITING FUNCTIONS **********/
void write_attr_text(hid_t loc_id,char *name,char * value){//write atribute text
	herr_t status;
	hid_t attr_id;
	const hsize_t dims[1] = {1};
	hid_t memtype;
	hid_t dataspace_id;
	char string[strlen(value)+10];
	strcpy(string,value);//I don't know why but I can not use the static string given in value for FIXED lenght strings, the enconding gets wrong, as I may not use H5T_VARIABLE I must copy to a variable string
	dataspace_id = H5Screate(H5S_SCALAR);
	
	memtype = H5Tcopy(H5T_C_S1);   
	H5Tset_strpad (memtype, H5T_STR_NULLTERM ); 
	H5Tset_cset(memtype,H5T_CSET_UTF8);
	H5Tset_size(memtype, strlen(string)+1 );// H5T_VARIABLE shall not be used here, section 3.1 prohibits this

	attr_id = H5Acreate(loc_id, name, memtype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attr_id, memtype, &string);
	status = H5Aclose(attr_id);
	status = H5Sclose(dataspace_id);
}

void write_attr_float(hid_t loc_id,char *name,float value){//write atribute float
	herr_t status;
	const hsize_t len=1;
	hid_t dataspace_id = H5Screate(H5S_SCALAR);
	hid_t attr_id=H5Acreate2(loc_id,name,H5T_IEEE_F32LE,dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attr_id,H5T_NATIVE_FLOAT, &value);
   	status = H5Aclose(attr_id);
	status = H5Sclose(dataspace_id); 
}
void write_attr_uint(hid_t loc_id,char *name,int value){//write atribute unsignet int
	herr_t status; 
	const hsize_t len=1;
	hid_t dataspace_id =  H5Screate(H5S_SCALAR);
	hid_t attr_id=H5Acreate2(loc_id,name,H5T_STD_U32LE,dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attr_id,H5T_NATIVE_UINT, &value);
   	status = H5Aclose(attr_id);
	status = H5Sclose(dataspace_id); 
}
void write_attr_double(hid_t loc_id,char *name,double value){//write atribute double
	herr_t status; 
	const hsize_t len=1;
	hid_t dataspace_id = H5Screate(H5S_SCALAR);
	hid_t attr_id=H5Acreate2(loc_id,name,H5T_IEEE_F64LE,dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attr_id,H5T_NATIVE_DOUBLE, &value);
   	status = H5Aclose(attr_id);
	status = H5Sclose(dataspace_id); 
}

void write_attr_long(hid_t loc_id,char *name,long int value){//write atribute long int 
	herr_t status; 
	const hsize_t len=1;
	hid_t dataspace_id = H5Screate(H5S_SCALAR);
	hid_t attr_id=H5Acreate2(loc_id,name,H5T_STD_I64LE,dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attr_id,H5T_NATIVE_LONG, &value);
   	status = H5Aclose(attr_id);
	status = H5Sclose(dataspace_id); 
}

/********* BUFFER MAKING FUNCTIONS **********/



unsigned char *make_buffer(int moment,Sweep * sweep,double max, double min){
	int nrays,nbins,i,j;
	nrays=sweep->h.nrays;
	nbins=sweep->ray[0]->h.nbins;
	float value;
	unsigned char scaled;
//	printf("%i  %f\n",moment,sweep->h.elev);
	unsigned char * buffer=(unsigned char *)malloc(sizeof(char)*nrays*nbins);
	if (buffer==NULL){printf("Memory Allocation Error. I'm leaving you"); exit(-1);}
	for(i=0;i<nrays;i++){
		for(j=0;j<nbins;j++){
			value=sweep->ray[i]->h.f(sweep->ray[i]->range[j]);
			if(value<131000){ //actually 131072 but, as one says, when it floats go safe
				scaled=(value-min)*(254)/(max-min)+1;
				if(scaled>255) scaled=255;
			}
			else
				scaled=0;
			buffer[i*nbins+j]=scaled;
		}
	}
	//printf("max %f, min %f\n",max,min);
	return buffer;

}
/********* ACCESSORY FUNCTIONS **********/

int get_unfolding(Ray * ray){
	if (ray->h.prf2<1) return 0;
	float frac=ray->h.prf/ray->h.prf2;
	if(frac>17/12) return 1;
	if(frac>31/24) return 2;
	if(frac>9/8) return 3;
	return 0;
}

unsigned long long int unix_time(Ray_header header){
	time_t time;
	struct tm tmTime;
	double time2;
	unsigned long long int time3;

	tmTime.tm_sec=0;	
	tmTime.tm_min=header.minute	;
	tmTime.tm_hour=header.hour	;
	tmTime.tm_mday=header.day	;
	tmTime.tm_mon=header.month-1;
	tmTime.tm_year=header.year-1900;	
	tmTime.tm_wday=-1	;
	tmTime.tm_yday=-1;
	tmTime.tm_isdst=-1;

	time=mktime(&tmTime);
	time=time-timezone;
	time2=(((double)time)+header.sec);
	time3=time2*pow(10,6);
	return time3;

}



double get_max(Sweep *sweep){

	int i,j;
	float max,current=sweep->ray[0]->h.f(sweep->ray[0]->range[0]);
	double ret;
	max=-131000;
	for(i=0;i<sweep->h.nrays;i++){
		for(j=0;j<sweep->ray[i]->h.nbins;j++){
			current=sweep->ray[i]->h.f(sweep->ray[i]->range[j]);
			if(max<current && current<131000) max=current;
		}
	}
	ret=max;
	return ret;
}

double get_min(Sweep *sweep){

	int i,j;
	float min,current=sweep->ray[0]->h.f(sweep->ray[0]->range[0]);
	double ret;
	min=131000;
	for(i=0;i<sweep->h.nrays;i++){
		for(j=0;j<sweep->ray[i]->h.nbins;j++){
			current=sweep->ray[i]->h.f(sweep->ray[i]->range[j]);
			//if(current<131000) printf("current %f  sweep %f ray %i range %i \n",current,sweep->h.elev,i,j);
			if(min>current) min=current;
		}
	}
	ret=min;
	return min;
}



void get_index_first_ray_in_sweep(Sweep *sweep,int *ray_index) // get the index of the ray with smallest time in sweep;
{
	int i,idx=-1,temp_int;
	unsigned long long int now,old=-1;
	for (i=0;i<sweep->h.nrays;i++){
		temp_int=sweep->ray[i]->h.sec;
		now=sweep->ray[i]->h.year*10000000000+sweep->ray[i]->h.month*100000000+sweep->ray[i]->h.day*1000000+sweep->ray[i]->h.hour*10000+sweep->ray[i]->h.minute*100+temp_int;
		if(now<old){old=now; idx=i;}
		
	}
	*ray_index=idx;
}

void get_index_last_ray_in_sweep(Sweep *sweep,int *ray_index)// get the index of the ray with smallest time in sweep;
{
	int i,idx,temp_int;
	unsigned long long int now,old=0;
	for (i=0;i<sweep->h.nrays;i++){
		temp_int=sweep->ray[i]->h.sec;
		now=sweep->ray[i]->h.year*10000000000+sweep->ray[i]->h.month*100000000+sweep->ray[i]->h.day*1000000+sweep->ray[i]->h.hour*10000+sweep->ray[i]->h.minute*100+temp_int;
		if(now>old){old=now; idx=i;}
	}
	*ray_index=idx;
} 
void get_index_first_ray_in_volume(Volume *volume,int *sweep_index,int *ray_index)// get the index of the ray with smallest time in volume;
{
	int i,sweep_idx,ray_idx,j,temp_int;
	unsigned long long int now,old=-1;
	for (j=0;j<volume->h.nsweeps;j++){
		get_index_first_ray_in_sweep(volume->sweep[j],&i);
		temp_int=volume->sweep[j]->ray[i]->h.sec;
		now=volume->sweep[j]->ray[i]->h.year*10000000000+volume->sweep[j]->ray[i]->h.month*100000000+volume->sweep[j]->ray[i]->h.day*1000000+volume->sweep[j]->ray[i]->h.hour*10000+volume->sweep[j]->ray[i]->h.minute*100+temp_int;
		if(now<old){old=now; sweep_idx=j;ray_idx=i;}
	}
	*ray_index=ray_idx;
	*sweep_index=sweep_idx;
}
void get_index_last_ray_in_volume(Volume *volume,int *sweep_index,int *ray_index)// get the index of the ray with smallest time in volume;
{
	int i,sweep_idx,ray_idx,j,temp_int;
	unsigned long long int now,old=0;
	for (j=0;j<volume->h.nsweeps;j++){
		get_index_last_ray_in_sweep(volume->sweep[j],&i);
		temp_int=volume->sweep[j]->ray[i]->h.sec;
		now=volume->sweep[j]->ray[i]->h.year*10000000000+volume->sweep[j]->ray[i]->h.month*100000000+volume->sweep[j]->ray[i]->h.day*1000000+volume->sweep[j]->ray[i]->h.hour*10000+volume->sweep[j]->ray[i]->h.minute*100+temp_int;
		if(now>old){old=now; sweep_idx=j;ray_idx=i;}
	}
	*ray_index=ray_idx;
	*sweep_index=sweep_idx;
}
