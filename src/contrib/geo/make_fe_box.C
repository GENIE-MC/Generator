{
  gSystem->Load("libGeom");

  TGeoManager * gm = new TGeoManager("world","");

  TGeoMaterial * fe56_mat = new TGeoMaterial("fe56", 56, 26, 7.87);
  TGeoMedium *   fe56_med = new TGeoMedium  ("fe56", 1, fe56_mat);

  TGeoVolume * top = gm->MakeBox("Top",fe56_med,20000,20000,20000);
  gm->SetTopVolume(top);
  gm->CloseGeometry();

  TFile f("./boxg.root","recreate");
  gm->Write();
  f.Close();
}

