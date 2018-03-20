#define	path	"/mnt/Data/xuansun/MPM_geant4_2012-2013"

{
  TChain *chain = new TChain("anaTree");

  for(int i = 0; i < 500; i++)
  {
    chain->Add(Form("%s/Beta_PolE/analyzed_%i.root", path, i));

    if(i% 100 == 0)
    {
      cout << "Adding file " << i+1 << "/2000 from Beta_PolE." << endl;
    }

  }

  for(int i = 0; i < 500; i++)
  {
    chain->Add(Form("%s/Beta_PolW/analyzed_%i.root", path, i));

    if(i% 100 == 0)
    {
      cout << "Adding file " << i+1 << "/2000 from Beta_PolW." << endl;
    }

  }

  cout << "Merging files...." << endl;
  chain->Merge("mpm_sim_betas_500EWFiles_2012-2013.root");

}



