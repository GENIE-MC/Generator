/*
rootgINukeVal.C
genie Intranuke Data Analysis
Author Andrew Seel March 2011
updated Juan Manfredi July 2012

This program simplifies the production of graphs of energy/angle v. cross section for both GENIE generated events and published cross section data. This ROOT script is an upgrade to gINukeVal (circa 2008) that changes the format file format in order to allow for the graphing of N data files.

The script is called with the location of the format file, the directory of published data, the directory of GENIE data files, and the directory in which graphs should be saved. Do not put trailing slashes.

For example (Assuming rootgINukeVal.C lives above the EX directory)

root 'rootgINukeVal.C("EX/f-file","EX/PubData","EX/GenieFiles","EX/PngFiles")'
*/


/*
The new format file is tag based. There are three tags.
  [RECORD]       <- starts a new graph
  [GENIE]        <- adds a new GENIE file to the graph
  [EXPERIMENTAL] <- adds a new published data file to the graph

The lines following each tag are fields which give instructions to the grapher on how to intepret the data, [RECORD] has 5, [GENIE] has 5, and [EXPERIMENTAL] only has 3 and one optional line.

FORMAT OVERVIEW:  (All white-space is ignored, each line is its own field, except for the xl,xu,yl,yu line, where each field is comma separated)

[RECORD]
  xl,xu,yl,yu
  Graph Type            <- Four valid values: Energy, Momentum, Angle, XS
  Graph Title
  Savename
  binFactor             <- Integer valued
[GENIE]
  filename
  cols
  cut
  dcth                  <- Ignored for Angle graphs
  legend title
[EXPERIMENTAL]
  filename
  cols
  legend title
  cut                   <- Normally uneccessary. Needed for XS filetypes
[GEANT]
  
*/


#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <TMath.h>
#include <TCanvas.h>
#include <TError.h>
#include <TTree.h>
#include <string>
#include <TCollection.h>
#include <TNtuple.h>
#include <vector>
#include <typeinfo>
#include <genieStyle.C>

//A small utility function that cuts off wrapping spaces
string trim(string in){
  size_t f = in.find_first_not_of(' ',0);
  size_t l = in.find_last_not_of(' ',0);
  if(l-f-1 > f){
    in = in.substr(f,l-f-1);
  }
  return in;
}

//Holds data associated with each .ginuke or .txt file
//  Is responsible for reading the file and generating the associated TTree
class DataFile
{
public:
  string dir;
  string filename;
  string name;
  string title;
  string cols;
  string gType;
  string cut;
  float dcth;
  int color;
  TFile* file;
  bool isValid; //=false
  bool GetData();
  DataFile();
  DataFile(string dname, string dtitle);
  DataFile(string dtype,string dir, string ddname, string dtitle,string dcols,string dcut,double ddcth, int dcolor);
  string getName(){return name;};
  int GetID() {return this->ID;};
  bool valid() {return isValid;};
  int analyzeTextFile();
  TTree* dataTree;
  TNtuple* dataTuple;
};

DataFile::DataFile(){}
DataFile::DataFile(string dname, string dtitle){
  name = dname;
  title = dtitle;
  isValid = GetData();
}

DataFile::DataFile(string dtype,string ddir, string ddname, string dtitle, string dcols, string dcut, double ddcth, int dcolor){
  dir = ddir;
  filename = ddname;
  gType = dtype;
  name = ddir+ddname;
  title = dtitle;
  cols = dcols;
  cut = dcut;
  dcth = ddcth;
  color = dcolor;
  if (color>=5) {color++;}
  isValid = GetData();
}

DataFile::GetData(){
  stringstream convert;
  int nameVal = 0;
  convert.str(filename);
  convert >> nameVal;
  bool success = false;
  //Find file type
  string type = "";
  size_t pos = filename.find_last_of('.');
  if(gType.compare("XS")==0){
    TNtuple newData ("newData","", cols.c_str() );
    newData.ReadFile(name.c_str());
    cout<<"TNtuple created"<<endl;
    dataTuple = (TNtuple*) newData.Clone();
    cout<<name.c_str()<<endl;
    cout<<dataTuple<<endl;
    success = true;
  }
  else{
    if(pos != string::npos){
      if(filename[pos+1] == 'r'){
        if(filename[pos-1] == 'e'){
          type.assign("ginuke");
        }
        else{
          type.assign("gst");
        }
      }
      else{
        type.assign("txt");
      }
    }
    else{
      type.assign("chain");
    }
    if(type[0]!='t'){
      if(type[0]=='g'){
        //Normal .ginuke file, Easy
        file = new TFile(name.c_str());
        if(file->IsOpen()){
          TTree* dummy = file->Get(type.c_str());
          dataTree = (TTree*)dummy->Clone(); 
          delete dummy;
          success=true;
        }
        else{
          cout<<name.c_str()<<" was not found. Skipping [GENIE] tag."<<endl;
          success=false;
        }
      }
      //Chaining together some files
      else{
        TChain chain("ginuke");
        TChain* dummy = &chain;
        //Need to have name separate from directory in order to properly construct this filename.
        string n = dir+"/gntp.*"+filename+"*.ginuke.root";
        cout<<"Making a chain called "<<n<<endl;
        chain.Add(n.c_str());
        cout<<"Chain created."<<endl;
        dataTree = (TTree*)dummy->Clone();
        cout<<"dummy cloned"<<endl;
        success=true;
      }
    }
    else{
      //Grabbing simple data from the text file
      dataTree = new TTree("textData",title.c_str());
      int nCol = this->analyzeTextFile();
      success = true;
      if(gType.compare("Energy")==0){
        if(nCol==3){dataTree->ReadFile(name.c_str(),"E/D:xsec:err1");}
        else if(nCol==2){dataTree->ReadFile(name.c_str(),"E/D:xsec");}
        else{cout<<"Published data file has an unrecognized format."<<endl;
             success=false;}
      }
      if(gType.compare("Momentum")==0){
        if(nCol==3){dataTree->ReadFile(name.c_str(),"ph/D:xsec:err1");}
        else if(nCol==2){dataTree->ReadFile(name.c_str(),"ph/D:xsec");}
        else{cout<<"Published data file has an unrecognized format."<<endl;
             success=false;}
      }
      if(gType.compare("Angle")==0){
        if(nCol==3){dataTree->ReadFile(name.c_str(),"cth/D:xsec:err1");}
        else if(nCol==2){dataTree->ReadFile(name.c_str(),"cth/D:xsec");}
        else{cout<<"Published data file has an unrecognized format."<<endl;
             success=false;}
      }
    }
  }

  return success;
}


int DataFile::analyzeTextFile(){
  //Determines the number of columns of numbers in the text file
  // Assumes columns are separated by spaces
  // Only actually determines if there are more than two columns, if so, it assumes there are 3
  int size=0;
  ifstream dataStream(name.c_str(),ios::in);
  string line;
  if(dataStream.is_open()){
    getline(dataStream,line);
    line = trim(line);
    while(line[0]=='#'){
      getline(dataStream,line);
    } 
    size_t begin = line.find(" ",0);
    size_t end = line.find_first_not_of(" ",begin);
    size_t again = line.find(" ",end);
    if(again!=string::npos){
      size=3;
    }
    else{
      size=2;
    }
  }
  dataStream.close();
  return size;
}

//Here resides all of the information contained in the current record.
//  It is responsible for parsing the format file for each record, and generating DataFile objects
//  when asked nicely.
class FormatFile{
public:
  string fFile;
  string xl,xu,yl,yu;
  vector<string> fileNames;
  vector<string> cols;
  vector<string> cuts;
  vector<double> dcths;
  vector<string> titles;
  vector<string> dataSource;
  string mtitle;
  string savename;
  string type;
  bool logx = false;
  bool logy = false;
  int binFactor = 1;
  FormatFile(string gtype);
  bool process(size_t record);
  int numFiles(){return fileNames.size();};
  DataFile* makeDataFile(int which,string dir=".");
  string fetchGTitle(){return mtitle;};
  int numOfType(string type);
};

FormatFile::FormatFile(string fileName){
  fFile = fileName;
}

int FormatFile::numOfType(string type){
  int count = 0;
  int i;
  for(i=0;i<dataSource.size();i++){
    if(dataSource[i].compare(type)==0){
      count++;
    }
  }
  return count;
}

bool FormatFile::process(size_t record){
  //Clearing variables to prevent contamination from previous process runs
  fileNames.clear();
  cols.clear();
  cuts.clear();
  dcths.clear();
  titles.clear();
  dataSource.clear();
  logx=false;
  logy=false;
  binFactor = 1;
  mtitle = "";
  savename = "";


  bool good = false;
  ifstream formatStream(fFile.c_str(),ios::in);
  string line;
  vector<string> tokens;
  string temp;
  if(formatStream.is_open()){
    //Get the recordth record tag 
    int depth = 0;
    size_t i=0;
    while(i<record && formatStream.eof()==false){
      getline(formatStream,line);
      line = trim(line);
      //Searching the current line for the tag markers
      size_t open = line.find_first_of('[',0);
      if(open!=string::npos){
        size_t close = line.find_first_of(']',0);
        if(close <= open){
          cout<<"Malformed Tag on following line"<<endl;
          cout<<"  "<<line<<endl;
          i = record+2;
        }
        else{
          //Pulling out tag contents
          string curTag = line.substr(open+1,close-open-1);
          if(curTag.compare("RECORD")==0){
            i++;
          }
        }
      }
    }
    //If there were enough records
    if(i == record){
      cout<<"Record "<<record<<" found, processing."<<endl;
      curTag = "PROCESS";
      good = true;
      depth = 1;
    }

    //If we hit another RECORD tag or the end of the file, we're done
    while(good == true && curTag.compare("RECORD")!=0 && formatStream.eof()==false){
      getline(formatStream,line);
      if(line[0]!='#'){ //Ignore lines starting with a hash
        line = trim(line);
        //Check if line is the opening of a new section, then extract the tag
        size_t open = line.find_first_of('[',0);
        if(open!=string::npos){
          size_t close = line.find_first_of(']',0);
          curTag = line.substr(open+1,close-open-1);
          depth = 0;
        }
        //Getting data out of the RECORD
        if(curTag.compare("PROCESS")==0){
          if(depth==1){
            //Get the coordinates of the graph
            size_t commaPos = line.find_first_of(',',0);
            xl = line.substr(0,commaPos);
            size_t nextPos = line.find_first_of(',',commaPos+1);
            xu = line.substr(commaPos+1,nextPos-commaPos-1);
            commaPos = nextPos;
            nextPos = line.find_first_of(',',commaPos+1);
            yl = line.substr(commaPos+1,nextPos-commaPos-1);
            commaPos = nextPos;
            nextPos = line.find_first_of(',',commaPos+1);
            yu = line.substr(commaPos+1,nextPos-commaPos-1);
          }
          if(depth==2){type = line;}
          if(depth==3){mtitle = line;}
          if(depth==4){savename = line;}
          if(depth==5){
            stringstream iss;
            iss.str(line); iss.clear();
            iss >> binFactor;
          }
        }
        else if(curTag.compare("GENIE")==0){ 
          if(depth==0){dataSource.push_back(curTag);}
          if(depth==1){fileNames.push_back(line);}
          if(depth==2){cols.push_back(line);}
          if(depth==3){cuts.push_back(line);}
          if(depth==4){
            stringstream iss;
            iss.str(line); iss.clear();
            double val;
            iss >> val;
            dcths.push_back(val);
          }
          if(depth==5){titles.push_back(line);}
        }
        else if(curTag.compare("EXPERIMENTAL")==0){
          if(depth==0){dataSource.push_back(curTag);}
          if(depth==1){fileNames.push_back(line);}
          if(depth==2){
            //Give some defaults to the stacks
            dcths.push_back(1);
            if(type.compare("XS")!=0){
              cuts.push_back("");
            }
            //Push a valid datum onto another stack
            cols.push_back(line);
          }
          if(depth==3){titles.push_back(line);}
          if(depth==4){cuts.push_back(line);}
        }
        else if(curTag.compare("RECORD")!=0){
          //Someone messed up the format file
          cout<<"Bad Tag: " <<line<<endl;
          good = false;
        }
        depth++;
      }
    }
    if(type.compare("Angle")!=0){
      cout<<"Nonangular distribution"<<endl;
      if(yu[0] == '-'){logy = true; yu = yu.substr(1);}
      if(yl[0] == '-'){logy = true; yl = yl.substr(1);}
      if(xu[0] == '-'){logx = true; xu = xu.substr(1);}
      if(xl[0] == '-'){logx = true; xl = xl.substr(1);}
    }
      //cout<<"("<<xl<<","<<yl<<"),("<<xu<<","<<yu<<")"<<endl;
  }
  return good;
}

DataFile* FormatFile::makeDataFile(int which,string dir,string source){
  int num=0;
  int i=0;
  //Find the position of the desired file in the vectors
  while(num <= which && i < dataSource.size()){
    if(dataSource[i].compare(source)==0){ num++; }
    i++;
  }
  i--;
  dir = dir+"/";
  //Setting some default values
  if(this->type.compare("Angle")==0){
    DataFile* temp = new DataFile(type,dir,fileNames[i], titles[i],cols[i],cuts[i],1.0,i+1);
  }
  else if(this->type.compare("Momentum")==0||this->type.compare("XS")==0){
    DataFile* temp = new DataFile(type,dir,fileNames[i], titles[i],cols[i],cuts[i],dcths[i],i+1);
  }
  else{
    DataFile* temp = new DataFile(type,dir,fileNames[i], titles[i],cols[i],cuts[i],dcths[i],i+1);
  }
  return temp;
}






//**************
//MAIN METHOD
//**************
//  Loops through the fFile, generating (and saving) a graph for each Record
//    Will look for EXPERIMENTAL files in dataDir, GENIE files in ROOTDir, and will
//    place all saved .png files into saveDir. The directories should not have trailing slashes
int rootgINukeVal(char* fFile, char* dataDir = ".", char* ROOTDir = ".",char* saveDir=".")
{
  string tFile (fFile);
  string dDir (dataDir);
  int legendSize; //counting entries in legend
  float TextSize, y1; //adjusted font size and lower bound on legend
  FormatFile format (tFile);
  TCanvas* cans;
  set_root_env();
  size_t curRecord=1;
  bool doCurrent = false;
  doCurrent = format.process(curRecord);
  if(doCurrent == false){
    //This if statement shouldn't be necessary, but if it
    // isn't here, the while loop tries to execute anyways...
   cout<<"Format File not found"<<endl;
   return 0;
  }
  while(doCurrent == true){
    cout<<format.numOfType("GENIE")<<" root files and "<<format.numOfType("EXPERIMENTAL")<<" published data files "<<endl;

    int legendSize = format.numOfType("GENIE") + format.numOfType("EXPERIMENTAL");
    stringstream convert;
    //gStyle->SetErrorX(0);
    //gStyle->SetStyle(genieStyle);
    set_root_env();
    Double_t size = 0;
    Int_t num = 0;
    Int_t target = 0;
    Int_t A = 0;
    Double_t factor = 0;
    TH1F* htemp = 0;
    bool good = false;
    //gROOT->SetStyle("T2K");

    //the FormatFile class holds the graph limits as strings, converting here 
    Double_t xl, yl, xu, yu;
    convert.str(format.xl); convert.clear();
    convert >> xl;
    convert.str(format.yl); convert.clear();
    convert >> yl;
    convert.str(format.xu); convert.clear();
    convert >> xu;
    convert.str(format.yu); convert.clear();
    convert >> yu;

    //Creating a new canvas and setting some parameters 
    string canName;
    canName.assign(curRecord,'*');
    cans= new TCanvas(canName.c_str(),format.fetchGTitle().c_str());
    cans->cd();
    TPad* curP = gPad;
    if(format.logx==true){gPad->SetLogx(1);}
    if(format.logy==true){gPad->SetLogy(1);}
    
    TextSize = .035;
    y1 = 1.05-(.075*legendSize);//1, .075

    TLegend* leg1 = new TLegend(.52,y1,1,1,"");//.6
    leg1->SetTextSize(TextSize);
    //TLegend* leg1 = new TLegend(.6,.8,1,1,"");

    //Get the frame, set parameters, redraw frame
    TH1F* hf1 = (TH1F*) cans->DrawFrame(xl,yl,xu,yu);
    hf1->SetTitle(format.mtitle.c_str());
    if(format.type.compare("Angle")==0){
      hf1->GetXaxis()->SetTitle("cos(#theta)");
      hf1->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} [#frac{mb}{sr}]");
    }
    else if(format.type.compare("Momentum")==0){
      hf1->GetXaxis()->SetTitle("Momentum [Mev]");
      hf1->GetYaxis()->SetTitle("#frac{d#sigma}{dp} [#frac{mb}{MeV}]");
    }
    else if(format.type.compare("XS")==0){
      hf1->GetXaxis()->SetTitle("Energy [MeV]");
      hf1->GetYaxis()->SetTitle("#sigma (mb)");
    }
    else{
      hf1->GetXaxis()->SetTitle("Energy [MeV]");
      hf1->GetYaxis()->SetTitle("#frac{d#sigma}{d#OmegadE} [#frac{mb}{sr#upointMev}]");
    }
    //hf1->GetXaxis()->SetNdivisions(-50202);
    hf1->GetYaxis()->CenterTitle();
    hf1->Draw();


    //Loop over each GINUKE file, scaling and drawing each
    int numRoots = format.numOfType("GENIE");
    int k;
    int markerStyle=20;

 
    for(k=0;k<numRoots;k++){
      DataFile* simData = format.makeDataFile(k,ROOTDir,"GENIE");
      if(format.type.compare("XS")==0){
        cans->cd();
        //Do things completely differently
        int curCol = simData->color;
        simData->dataTuple->SetMarkerColor(curCol);
        simData->dataTuple->SetMarkerStyle(markerStyle);
	markerStyle++;
        simData->dataTuple->SetLineStyle(2);//2
        simData->dataTuple->SetLineColor(curCol);
        cout<<"About to draw tuple"<<endl;
        //simData->dataTuple->Draw();
        cout<<simData->cut.c_str()<<endl;
        simData->dataTuple->Draw(simData->cut.c_str(),"","line psame L");
	leg1->AddEntry(simData->dataTuple,simData->title.c_str(),"P");
        cout<<"Tuple drawn"<<endl;
      }
      else{
        TCanvas* tempVas = new TCanvas("tempName","No title");
        tempVas->cd();
        //TPad* curP = gPad;
        string newCut = simData->cut;
        newCut = newCut +"&&"+simData->cols+"<="+format.xu;
        newCut = newCut + "&&"+simData->cols+">"+format.xl;
        simData->dataTree->Draw(simData->cols.c_str(),newCut.c_str(), "L");
        //tempVas is used in order to not clobber cans's htemp
        if(simData->valid()){

          //Grab the associated histogram from tempVas, apply scaling factor
          htemp = (TH1F*) gPad->GetPrimitive("htemp");
          TH1F* hist1;
	  TH1F* hist1error;
          good=false;
          //Make sure that something actually existed in the cut
          if (htemp != 0x0) {
            good=true;
            hist1 = (TH1F*) htemp->Clone("hist1");
            size = hist1->GetBinWidth(1);
            hist1->Sumw2();
            hist1->Rebin(format.binFactor);
            num = simData->dataTree->GetEntries();
            simData->dataTree->SetBranchAddress("tgt",&target);
            simData->dataTree->GetEntry(1);
            A = (target/10) - (target/10000)*1000;
            factor = TMath::Power(3*1.4*TMath::Power(A,(1.0/3.0)),2)*10.0/(2*num*simData->dcth);
            hist1->Scale(factor,"width");
            int curCol = simData->color;
            hist1->SetMarkerColor(curCol);
            hist1->SetMarkerStyle(markerStyle);
	    markerStyle++;
            hist1->SetLineStyle(2);
            hist1->SetLineColor(curCol);
          }
          else{
            cout<<"Nothing was found in the cut of "<<simData->filename<<endl;
          }
          if(good){
            //Draw a copy of histogram to cans on top of any other histograms already there
            cans->cd();
            //TPad* curP = gPad;
	    leg1->AddEntry(hist1,simData->title.c_str());
            hist1->DrawCopy("e1 psame");//e1 psame
	    hist1->Draw("hist l same");
	    
	    //make and draw an extra histogram to have solid error bars
	    hist1error = (TH1F*) hist1->Clone("hist1error"); 
	    hist1error->SetLineStyle(1);
	    hist1error->DrawCopy("e1 same");
          }
        }
        else{
          cout<<"Something is wrong with data file "<<simData->filename<<endl;
          cans->cd();
        }
        tempVas->Close();
      }
    }

    //Draw all of the Experimental files, I'm pretty sure that this will only handle
    // the simple tree type generated from the .txt type data files.
    cans->cd();
    TPad* curP = gPad;
    int numFiles = format.numOfType("EXPERIMENTAL");

    /*
    int j = 0;
    
    for(j=0;j<3;j++){
      DataFile* experimental = format.makeDataFile(0,dDir,"EXPERIMENTAL");
      TGraphErrors* data1;
      if(format.type.compare("XS")==0){
        experimental->dataTuple->Draw(experimental->cut.c_str(),"","goff");
        data1 = new TGraphErrors(experimental->dataTuple->GetSelectedRows(),experimental->dataTuple->GetV1(), experimental->dataTuple->GetV2(),experimental->dataTuple->GetV3(),experimental->dataTuple->GetV4());
      }
      else{
        experimental->dataTree->Draw(experimental->cols.c_str(),"","goff");
        data1 = new TGraphErrors(experimental->dataTree->GetSelectedRows(),experimental->dataTree->GetV2(), experimental->dataTree->GetV1(),0,experimental->dataTree->GetV3());
      }
      //data1->SetLineStyle(0);
      data1->SetMarkerColor(experimental->color);
      data1->SetMarkerStyle(markerStyle);
      markerStyle++;
      leg1->AddEntry(data1,experimental->title.c_str(),"P");
      data1->Draw("p same");
    }
    */
    
    for(j=0;j<numFiles;j++){
      DataFile* experimental = format.makeDataFile(j,dDir,"EXPERIMENTAL");
      TGraphErrors* data1;
      if(format.type.compare("XS")==0){
        experimental->dataTuple->Draw(experimental->cut.c_str(),"","goff");
        data1 = new TGraphErrors(experimental->dataTuple->GetSelectedRows(),experimental->dataTuple->GetV1(), experimental->dataTuple->GetV2(),experimental->dataTuple->GetV3(),experimental->dataTuple->GetV4());
      }
      else{
        experimental->dataTree->Draw(experimental->cols.c_str(),"","goff");
        data1 = new TGraphErrors(experimental->dataTree->GetSelectedRows(),experimental->dataTree->GetV2(), experimental->dataTree->GetV1(),0,experimental->dataTree->GetV3());
      }
      data1->SetLineStyle(3);
      data1->SetMarkerColor(experimental->color);
      data1->SetMarkerStyle(markerStyle);
      markerStyle++;
      leg1->AddEntry(data1,experimental->title.c_str(),"P");
      data1->Draw("p same");
    }
    
    
    //Draw the legend    

    //set_root_env();
    leg1->SetLineWidth(2);
    leg1->Draw();

    //add_plot_label("THIS IS A TEST #alpha #Alpha#gamma", .5, .5, .05);

    //Save the record
    string saveName(saveDir);
    saveName = saveName+"/"+format.savename+".png";
    cans->SaveAs(saveName.c_str());

    //Try to find the next record
    curRecord++;
    cout<<"\nLooking for next record:"<<endl;
    doCurrent = format.process(curRecord);
  }
  cout<<"No more records found, exiting"<<endl;
  //delete [] entryArray;
  return 0;
}
