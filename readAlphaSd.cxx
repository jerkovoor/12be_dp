	///////////////////////////////////////
	/////////////S1 DETECTORS//////////////
	///////////////////////////////////////
	
	
	//find peaks and for the Sdr detectors
	double Sdrxp[1];
	//double Sdra[24],Sdrb[24],Sdrc[24];
	double Sdra[24];
	
	for ( int j=0; j<24; j++ ){
		nfound = s->Search ( hSd1r[j],2,"",0.1 ); //search for peaks in the hSdr detectors change the name of the histograme according to which peaks you want to search for Sd1r, Sd2r or Sur
       
		for ( int p=0; p<nfound; p++ ){Sdrxp[p] = txp[p];} 
        
		/*Sdra[j] = min ( min ( Sdrxp[0],Sdrxp[1] ), Sdrxp[2] );
		Sdrc[j] = max ( max ( Sdrxp[0],Sdrxp[1] ), Sdrxp[2] );

		if(Sdrxp[0]!=Sdra[j] && Sdrxp[0]!=Sdrc[j]) {Sdrb[j] = Sdrxp[0];} 
		else if(Sdrxp[1]!=Sdra[j] && Sdrxp[1]!=Sdrc[j]) {Sdrb[j] = Sdrxp[1];}
		else if(Sdrxp[2]!=Sdra[j] && Sdrxp[2]!=Sdrc[j]) {Sdrb[j] = Sdrxp[2];}  */
    }
  
	//open a txt file to store the positions of the peaks
	ofstream Sdr_alpha;
	Sdr_alpha.open ( "/home/jerome/12Be_exp/notarget_peak_position/Sd1r_notarget_peak.txt" ); //change the path and the name of the file acordingly
	//Sdr_alpha << " Strip / Peak 1 / 2 / 3  \n";
	Sdr_alpha << " Strip / Peak 1 \n";
	
	ofstream Sdr_alphaCounts;
	Sdr_alphaCounts.open ( "/home/jerome/12Be_exp/notarget_peak_position/Sd1r_notarget_counts.txt" ); //change the path and the name of the file acordingly
	//Sdr_alphaCounts << " Strip / Amplitude 1 / 2 / 3  \n";
	Sdr_alphaCounts << " Strip / Amplitude 1 \n";
	
	for ( int i=0; i<24; i++ ){ //loop over all the rings of S3 detectors to set the parameters and fit
		
		//if(i<6){
			/*fit_func3->SetParLimits ( 1,Sdra[i]-5,Sdra[i]+5 ); //remove here
			fit_func3->SetParLimits ( 4,Sdrb[i]-5,Sdrb[i]+5 );
			fit_func3->SetParLimits ( 6,Sdrc[i]-5,Sdrc[i]+5 );
			fit_func3->SetParLimits ( 2,0,10 );*/			//remove here
			
			fit_func1->SetParLimits ( 1,Sdrxp[i]-200,Sdrxp[i]+200 ); //2000,2500
			fit_func1->SetParLimits ( 2,0,200 ); //0,200
	//	}else if(i>5){
			//fit_func3->SetParLimits ( 1,Sdra[i]-3,Sdra[i]+3 );
			//fit_func3->SetParLimits ( 4,Sdrb[i]-3,Sdrb[i]+3 );
			//fit_func3->SetParLimits ( 6,Sdrc[i]-3,Sdrc[i]+3 );
			//fit_func3->SetParLimits ( 2,0,3 );
		//}
		hSd1r[i]->Fit ( "fit_func1","","",Sdrxp[i]-500,Sdrxp[i]+500 ); //change the histograme accordingly // 1000,3000

		//fill in the .txt files
		//Sdr_alpha << i << " "  << fit_func3->GetParameter ( 1 ) << " " << fit_func3->GetParameter ( 4 ) << " " << fit_func3->GetParameter ( 6 ) << endl;
		Sdr_alpha << i << " "  << fit_func1->GetParameter ( 1 ) << endl;
		// Sdr_alphaCounts << i << " " << fit_func3->GetParameter ( 0 ) << " " <<  fit_func3->GetParameter ( 3 ) << " " <<  fit_func3->GetParameter ( 5 ) << endl;
		Sdr_alphaCounts << i << " " << fit_func1->GetParameter ( 0 ) << endl;
	}//end of fit on Sd1r detector
    
    
  
	//find peaks and for the Sds detectors
	double Sdsxp[1];
	//double Sdsa[32],Sdsb[32],Sdsc[32],Sdsped[32];
	
	for ( int j=0; j<32; j++ ){ //loop over the sectors of the S3 detectors
		nfound = s->Search ( hSd1s[j],2,"",0.1 ); //change the name of the hitogram according to which detector you want to look at: Sd1s, Sd2s or Sus
		
		for ( int p=0; p<nfound; p++ ){Sdsxp[p] = txp[p];}
		
		/* Sdsa[j] = min ( min ( Sdsxp[0],Sdsxp[1] ), Sdsxp[2] );
		Sdsc[j] = max ( max ( Sdsxp[0],Sdsxp[1] ), Sdsxp[2] );

		if(Sdsxp[0]!=Sdsa[j] && Sdsxp[0]!=Sdsc[j]) {Sdsb[j] = Sdsxp[0];}
		else if(Sdsxp[1]!=Sdsa[j] && Sdsxp[1]!=Sdsc[j]) {Sdsb[j] = Sdsxp[1];}
		else if(Sdsxp[2]!=Sdsa[j] && Sdsxp[2]!=Sdsc[j]) {Sdsb[j] = Sdsxp[2];}     */  
                           
	}
  
	//open a txt file to store the positions of the peaks
	ofstream Sds_alpha;
	Sds_alpha.open ( "/home/jerome/12Be_exp/notarget_peak_position/Sd1s_notarget_peak.txt" );//change the path and the name of the file acordingly
	//Sds_alpha << " Strip / Peak 1 / 2 / 3  \n";
	Sds_alpha << " Strip / Peak 1 \n";
  
	ofstream Sds_alphaCounts;
	Sds_alphaCounts.open ( "/home/jerome/12Be_exp/notarget_peak_position/Sd1s_notarget_counts.txt" );//change the path and the name of the file acordingly
	//Sds_alphaCounts << " Strip / Amplitude 1 / 2 / 3  \n";
	Sds_alphaCounts <<  " Strip / Amplitude 1 \n";
    
	for ( int i=0; i<32; i++ ){ //loop over the sectors of S3 detectors to see the parameters and fit
		//if(i<16){
			/*fit_func3->SetParLimits ( 1,Sdsa[i]-5,Sdsa[i]+5 );
			fit_func3->SetParLimits ( 4,Sdsb[i]-5,Sdsb[i]+5 );
			fit_func3->SetParLimits ( 6,Sdsc[i]-5,Sdsc[i]+5 );
			fit_func3->SetParLimits ( 2,0,10 );*/
			fit_func2->SetParLimits ( 1,Sdsxp[i]-200,Sdsxp[i]+200 );
			fit_func2->SetParLimits ( 2,0,200 ); 
		//}else if(i>15){
			//fit_func3->SetParLimits ( 1,Sdsa[i]-3,Sdsa[i]+3 );
			//fit_func3->SetParLimits ( 4,Sdsb[i]-3,Sdsb[i]+3 );
			//fit_func3->SetParLimits ( 6,Sdsc[i]-3,Sdsc[i]+3 );
			//fit_func3->SetParLimits ( 2,0,3 );
		//}
  
		hSd1s[i]->Fit ( "fit_func2","","",Sdsxp[i]-500,Sdsxp[i]+500 ); //change the histograme name acordingly

		//fill in the .txt files
		//Sds_alpha << i << " "  << fit_func3->GetParameter ( 1 ) << " " << fit_func3->GetParameter ( 4 ) << " " << fit_func3->GetParameter ( 6 ) << endl;
		Sds_alpha << i << " "  << fit_func2->GetParameter ( 1 ) << endl;
		//Sds_alphaCounts << i << " " << fit_func3->GetParameter ( 0 ) << " " <<  fit_func3->GetParameter ( 3 ) << " " <<  fit_func3->GetParameter ( 5 ) << endl;
		Sds_alphaCounts << i << " " << fit_func2->GetParameter ( 0 ) <<  endl;
	}//end of fit on Sd1s detector