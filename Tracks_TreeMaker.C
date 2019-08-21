
gROOT->ProcessLine(".L common/gallery_includes.h+O");
gROOT->ProcessLine(".L common/larsoftobj_includes.h+O");
gROOT->ProcessLine(".L common/My_ROOT_Utilities.h");

//setup our input file
std::vector<std::string> filenames;
std::string file_path;
size_t max_files;

//setup our output file
TFile *f_output;

//setup for our input tags
art::InputTag tag_acpt_intime;
art::InputTag tag_tracks;
art::InputTag tag_calo;
art::InputTag tag_hits;
art::InputTag tag_hit_truth_assn;

//filenames = GetFileNameList("/uboone/app/users/jaz8600/work/MCC9_1/DetectorSystematics/data_extbnb_mcc9.1_v08_00_00_16_18_run1_reco2.list","");
filenames = GetFileNameList("/uboone/app/users/jaz8600/work/MCC9_1/DetectorSystematics/AnodeCrossersIso.list","");

//clip the filename list
max_files = 10;
if(filenames.size()>max_files)
    filenames.resize(max_files)

out = new TFile("outfile_MC.root","RECREATE");

tag_acpt_intime = { "acpttrigtagger" };
tag_tracks = { "pandora" };
tag_calo = { "pandoracaliSCE" };
tag_hits = { "gaushit" };
tag_hit_truth_assn = { "gaushitTruthMatch" };

TTree* fTree;
int Run;
int Subrun;
int Event;

int Ntrks;
int Nacpt;
int n_acpt;
int Nhits;
int n_hits;

double hit_Q;
double hit_A;
double hit_sigma;
double hit_time;

int channel;
double wire_Q;
double wire_A;
double wire_sigma;
double wire_time;
int wire_begin;
int wire_end;

int hit_plane;
int hit_isMC;

double trk_x;
double trk_y;
double trk_z;
double trk_L;
double trk_cosTheta;
double trk_phi;
double trk_ThetaXZ;
double trk_ThetaYZ;
double trk_dEdx;
double trk_fracMC;

fTree = new TTree("hit_v_track","HitPropertiesTrackProperties");

fTree->Branch("Run",&Run);
fTree->Branch("Subrun",&Subrun);
fTree->Branch("Event",&Event);
fTree->Branch("Ntrks",&Ntrks);
fTree->Branch("Nacpt",&Nacpt);  
fTree->Branch("n_acpt",&n_acpt);  
fTree->Branch("Nhits",&Nhits);  
fTree->Branch("n_hits",&n_hits);  
fTree->Branch("hit_Q",&hit_Q);  
fTree->Branch("hit_A",&hit_A);  
fTree->Branch("hit_sigma",&hit_sigma);  
fTree->Branch("hit_time",&hit_time);  
fTree->Branch("hit_plane",&hit_plane);  
fTree->Branch("hit_isMC",&hit_isMC);

fTree->Branch("channel",&channel);
fTree->Branch("wire_Q",&wire_Q);
fTree->Branch("wire_A",&wire_A);
fTree->Branch("wire_sigma",&wire_sigma);
fTree->Branch("wire_time",&wire_time);
fTree->Branch("wire_begin",&wire_begin);
fTree->Branch("wire_end",&wire_end);

fTree->Branch("trk_x",&trk_x);  
fTree->Branch("trk_y",&trk_y);  
fTree->Branch("trk_z",&trk_z);  
fTree->Branch("trk_L",&trk_L);  
fTree->Branch("trk_cosTheta",&trk_cosTheta);  
fTree->Branch("trk_phi",&trk_phi);  
fTree->Branch("trk_ThetaXZ",&trk_ThetaXZ);  
fTree->Branch("trk_ThetaYZ",&trk_ThetaYZ);  
fTree->Branch("trk_dEdx",&trk_dEdx);  
fTree->Branch("trk_fracMC",&trk_fracMC); 

fTree->Clear();

TF1 *fgaus_norm,*fgaus_mod;

fgaus_norm = new TF1("fgaus_norm","gaus",-20,20);
fgaus_norm->SetParameters(1.0,0.0,1.0);
fgaus_mod = new TF1("fgaus_mod","gaus",-20,20);
fgaus_mod->SetParameters(1.0,0.0,1.0);

%%cpp -d
void CalcWvfmCenterAndRMS(std::vector<float> const& in_wvfm,
                         float & center,
                         float & rms,
                         float & integral)
{
    rms=0;
    center = 0;
    integral = 0;
    for(size_t i_t=0; i_t<in_wvfm.size(); ++i_t){
        center += in_wvfm[i_t]*i_t;
        integral += in_wvfm[i_t];
    }
    
    center = center / integral;
    for(size_t i_t=0; i_t<in_wvfm.size(); ++i_t)
        rms += in_wvfm[i_t]*(i_t-center)*(i_t-center);
    rms = std::sqrt(rms/integral);
    
}

/*
                                wire_begin = range.begin_index();
                                wire_end = wire_begin+range.size();
                                wire_Q = 0;
                                wire_A = 0;
                                wire_time = 0;
                                wire_sigma = 0;
                                for(size_t i_t = 0; i_t<range.data().size(); ++i_t){
                                    auto val = range.data().at(i_t);
                                    auto tick = wire_begin+i_t;
                                    if(val>wire_A) wire_A = val;
                                    wire_Q += val;
                                    wire_time += val*tick;
                                }
                                wire_time = wire_time/wire_Q;
                                for(size_t i_t = 0; i_t<range.data().size(); ++i_t){
                                    auto val = range.data().at(i_t);
                                    auto tick = wire_begin+i_t;
                                    wire_sigma += val*(tick-wire_time)*(tick-wire_time);
                                }
                                wire_sigma = std::sqrt(wire_sigma/wire_Q);
*/


%%cpp -d

void SmearWaveform(std::vector<float> const& in_wvfm, std::vector<float> & out_wvfm,
                   float f_s, unsigned int mask_extent=2){

    fgaus_norm->SetRange(0,in_wvfm.size());
    fgaus_mod->SetRange(0,in_wvfm.size());

    float in_center,in_sigma,in_integral;
    
    CalcWvfmCenterAndRMS(in_wvfm,in_center,in_sigma,in_integral);
    
    fgaus_norm->SetParameter(1,in_center);
    fgaus_norm->SetParameter(2,in_sigma);
    float integral_norm = fgaus_norm->Integral(0,in_wvfm.size());
    
    fgaus_mod->SetParameter(1,in_center);
    fgaus_mod->SetParameter(2,in_sigma*f_s);
    float integral_mod = fgaus_mod->Integral(0,in_wvfm.size());

    /*    
    //create mask
    std::vector<float> mask(1+2*mask_extent);
    fgaus_norm->SetParameter(2,sigma_orig);
    fgaus_mod->SetParameter(2,f_s*sigma_orig);
    float mod_integral = fgaus_mod->Integral(-20,20);
    float norm_integral = fgaus_norm->Integral(-20,20);
    for(size_t i_m=0; i_m<mask.size(); ++i_m){
        float val = (float)i_m - float(mask_extent);
        mask[i_m] = (fgaus_mod->Eval(val)/mod_integral) / (fgaus_norm->Eval(val)/norm_integral);
        std::cout << " " << mask[i_m];
    }
    std::cout << std::endl;
    */
    out_wvfm.resize(in_wvfm.size());
    
    for(size_t i_t=0; i_t < in_wvfm.size(); ++i_t){
        std::cout << " " << fgaus_mod->Eval(i_t)/integral_mod << "/" << fgaus_norm->Eval(i_t)/integral_norm
                  << " (" << (fgaus_mod->Eval(i_t)/integral_mod) / (fgaus_norm->Eval(i_t)/integral_norm) << ")"
                  << std::endl;
        out_wvfm[i_t] = in_wvfm[i_t] * (fgaus_mod->Eval(i_t)/integral_mod) / (fgaus_norm->Eval(i_t)/integral_norm);
    }
    
}

std::vector<float> mywvfm;
std::vector<float> out_wvfm;

float my_ped,my_rms,my_integral;

//mywvfm = std::vector<float>({0,0,2,5,9,5,2,0,0});
/*
mywvfm = std::vector<float>(30);

fgaus_norm->SetParameter(1,15);
fgaus_norm->SetParameter(2,3);
fgaus_norm->SetParameter(0,10.0);

for(size_t i_t=0; i_t<mywvfm.size(); ++i_t)
{
    mywvfm[i_t] = fgaus_norm->Eval(i_t);
}

*/
mywvfm = std::vector<float>({0,0,0,0,10,10,10,10,0,0,0,0});


SmearWaveform(mywvfm,out_wvfm,0.5)

CalcWvfmCenterAndRMS(mywvfm,my_ped,my_rms,my_integral);
std::cout << "Original: " << my_ped << "," << my_rms << "," << std::endl;
CalcWvfmCenterAndRMS(out_wvfm,my_ped,my_rms,my_integral);
std::cout << "Modified: " << my_ped << "," << my_rms << "," << std::endl;

for (gallery::Event ev(filenames) ; !ev.atEnd(); ev.next()) {
        
    /// Prep our Branches
    Run = 0; //
    Subrun = 0; //
    Event = 0; //
    Ntrks = 0; //
    Nacpt = 0;//
    n_acpt = -1;//
    Nhits = 0;//
    n_hits = -1;//
    hit_Q = 0; //
    hit_A = 0; //
    hit_sigma = 0; //
    hit_time = 0; //
    hit_plane = 0; //
    trk_x = 0;//
    trk_y = 0;//
    trk_z = 0;//
    trk_cosTheta = 0;//
    trk_phi = 0;//
    trk_ThetaXZ = 0;//
    trk_ThetaYZ = 0;//
    trk_dEdx = 0;//
    trk_L = 0;//
    trk_fracMC = 0;
    hit_isMC = 0;
    
    channel = -1;
    wire_Q = 0;
    wire_A = 0;
    wire_sigma = 0;
    wire_time = 0;
    wire_begin = -1;
    wire_end = -1;
    
    //set run/subrun/event info
    Run = ev.eventAuxiliary().run();
    Subrun = ev.eventAuxiliary().subRun();
    Event = ev.eventAuxiliary().event();
    
    //We will now found the events that have a ACPT in-time track
    auto const& t0s = *(ev.getValidHandle<std::vector<anab::T0>>(tag_acpt_intime));
    
    //Skipping those that don't
    if(t0s.size() == 0) continue;
    Nacpt = t0s.size();
    
    //This associates the T0 tag with the track
    auto const &t0_assoc_handle =
      ev.getValidHandle<art::Assns<anab::T0, recob::Track>>(tag_acpt_intime);
    
    //Make a vector that will hold our tagged tracks
    std::vector<recob::Track> ACPT_tracks;

    // store the tagged tracks into that vector
    for(auto &assn : *t0_assoc_handle)
      ACPT_tracks.emplace_back(*(assn.second));
        
    // Now we'll need to set things up to collect the calorimetry data
    // Start by saying which tracks we want:
    auto const & track_list_ptr = ev.getValidHandle<std::vector<recob::Track> >(tag_tracks);
    auto const & track_list = (*track_list_ptr);
    
    Ntrks = track_list.size();
    
    //Next we'll get the associated calorimetries things:
    art::FindMany<anab::Calorimetry>  fmcal(track_list_ptr, ev, tag_calo);
    
    //This let's us find which hits are assocaited to a give track trajectory point
    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(track_list_ptr, ev, tag_tracks);
    
    //get hit assn to recob::Wire
    auto const& hit_handle = ev.getValidHandle<std::vector<recob::Hit>>(tag_hits);
    art::FindMany<recob::Wire> wires_per_hit(hit_handle,ev,tag_hits);

    //Find all the hits that are matched 
    // to MCParticles
    art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hit_handle,ev,tag_hit_truth_assn);
    bool is_simulation = particles_per_hit.isValid();    
    
    //Let's loop through our tracks and find our calorimetry things
    for(auto &trk : ACPT_tracks){      
        n_acpt++;
        for(int itrk = 0; itrk < int(track_list.size()); itrk++){
            if(trk.ID() == track_list.at(itrk).ID()){
                trk_L = trk.Length();
                
                // Now we have a track which is matched to a ACPT crossing track
                // now we want the hits for this track

                // This is the vector of hits
                auto vhit = fmthm.at(itrk);
                Nhits = vhit.size();

                if(is_simulation){
                
                    // Check fraction that are truth:
                    int n_truthhits = 0;
                    for(int ht = 0; ht < Nhits; ht++){
                        auto part_vec = particles_per_hit.at(vhit[ht].key());
                        if(part_vec.size() > 0)
                            n_truthhits++;
                    }

                    if(Nhits > 0)
                        trk_fracMC = double(n_truthhits)/double(Nhits);
                    else
                        trk_fracMC = 0.;
                }
                
                // This is the vector of traj point info
                // the Index() of this is the traj point of this track
                auto vmeta = fmthm.data(itrk);

                // Now we can get the calorimetry points
                std::vector<const anab::Calorimetry*> calos = fmcal.at(itrk);
                
                // this will count which calo point we're on
                int count = 0;
                
                //iterate through the planes:
                for(int pl = 0; pl < 3; pl++){

                    //iterate through track meta points :           
                    for(int vp = 0; vp < vmeta.size(); vp++){
                        // store the track trajectory point index
                        // for the track meta point
                        int ind = vmeta[vp]->Index();
                        
                        // check that the traj point is in the calorimetry point
                        // and belongs to the plane we are interested in 
                        if(track_list.at(itrk).HasValidPoint(ind) && vhit[vp]->WireID().Plane == pl){
                        
                            n_hits++;

                            if(is_simulation){
                                auto part_vec = particles_per_hit.at(vhit[vp].key());           
                                if(part_vec.size() > 0)
                                    hit_isMC = 1;
                            }
                            
                            // Grab the track traj point
                            // WE DON'T CURRENTLY USE THIS
                            // I kept it for testing purposes 
                            auto trjp = track_list.at(itrk).TrajectoryPoint(ind);
                
                            // Grab the calo point
                            auto calp = calos[pl]->XYZ()[count];
                            auto caldEdx = calos[pl]->dEdx()[count];
                
                            // We need to calculate the angles 
                            // of the calo points 
                            double Phi = 0;
                            double cosTheta = 0;
                            double ThetaXZ = 0;
                            double ThetaYZ = 0;
                            
                            auto angle = (count < vmeta.size()-1) ? 
                                (calos[pl]->XYZ()[count]) - (calos[pl]->XYZ()[count+1]) : 
                                (calos[pl]->XYZ()[count-1]) - (calos[pl]->XYZ()[count]);
                            
                            Phi = angle.Phi(); 
                            cosTheta = cos(angle.Theta());
                            ThetaXZ = atan2(angle.X(),angle.Z());
                            ThetaYZ = atan2(angle.Y(),angle.Z());           

                            // Grab the matched hit 
                            auto hit = vhit[vp];
                            
                            hit_A = hit->PeakAmplitude();
                            hit_Q = hit->Integral();
                            hit_sigma = hit->RMS();
                            hit_time = hit->PeakTime();
                            hit_plane = pl;
                
                            trk_x = calp.X(); 
                            trk_y = calp.Y(); 
                            trk_z = calp.Z(); 
                            trk_phi = Phi;
                            trk_cosTheta = cosTheta;
                            trk_dEdx = caldEdx;
                            trk_ThetaXZ = ThetaXZ;
                            trk_ThetaYZ = ThetaYZ;

                            channel = hit->Channel();
                            
                            //Grab wire
                            auto wire_vec = wires_per_hit.at(hit.key());

                            //there should be only one associate recob::Wire...
                            if(wire_vec.size()==0){
                                std::cout << "ERROR: No associated wire!" << std::endl;
                                continue;
                            }
                            else if(wire_vec.size()>1){
                                std::cout << "WARNING: More than one associated Wire? Only taking the first." << std::endl;
                            }
                            
                            auto const* wire_ptr = wire_vec[0];
                            auto const& signal_rois = wire_ptr->SignalROI();
                            for (auto const& range: signal_rois.get_ranges()) {
                                if(hit->PeakTime()<range.begin_index() || 
                                    hit->PeakTime()>range.begin_index()+range.size()) continue;
                                wire_begin = range.begin_index();
                                wire_end = wire_begin+range.size();
                                wire_Q = 0;
                                wire_A = 0;
                                wire_time = 0;
                                wire_sigma = 0;
                                for(size_t i_t = 0; i_t<range.data().size(); ++i_t){
                                    auto val = range.data().at(i_t);
                                    auto tick = wire_begin+i_t;
                                    if(val>wire_A) wire_A = val;
                                    wire_Q += val;
                                    wire_time += val*tick;
                                }
                                wire_time = wire_time/wire_Q;
                                for(size_t i_t = 0; i_t<range.data().size(); ++i_t){
                                    auto val = range.data().at(i_t);
                                    auto tick = wire_begin+i_t;
                                    wire_sigma += val*(tick-wire_time)*(tick-wire_time);
                                }
                                wire_sigma = std::sqrt(wire_sigma/wire_Q);
                                
                            }
                            /*
                            for(auto iROI = wire_ptr->begin_range(); iROI != wire_ptr->end_range(); ++iROI) {
                                auto const& ROI = *iROI;
                                //const int FirstTick = ROI.begin_index();
                                //const int EndTick = ROI.end_index();
                                //const float FirstADC = ROI[FirstTick];
                                
                                std::cout << ROI.size() << std::endl;
                                
                                //for (float ADC: ROI){}
                            } // for
                            */
                            fTree->Fill();
                                            
                            // this tracks the correct calorimetry 
                            // we are supposed to be anlayzing 
                            count++;
                            
                        }//end if desired calo point
                        
                    }//end loop over track metapoints
                    
                    count=0;
                    
                }//end loop over planes
                 
            }//end if the track we want
            
        }//end loop over all tracks
        
    }//end loop over ACPT_tracks
    
}//end loop over events

%jsroot on

TCanvas mycanvas("mycanvas","MyCanvas");

fTree->Draw("wire_time-hit_time");
mycanvas.Draw();

fTree->Draw("wire_Q-hit_Q","(wire_end-wire_begin)<50");
mycanvas.Draw();

fTree->Draw("wire_sigma");
mycanvas.Draw();

fTree->Draw("wire_sigma-hit_sigma");
mycanvas.Draw();
