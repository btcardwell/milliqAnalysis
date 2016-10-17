const char * dirName = "/data/users/MilliQan/caen/PeriodicDarkRate/";
const char * ext     = ".root";

void readAllFiles() {
    TSystemDirectory dir(dirName, dirName);
    TList * files = dir.GetListOfFiles();

    if (files) {
        TSystemFile * file;
        TString fName;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fName = file->GetName();
            if (!file->IsDirectory() && fName.EndsWith(ext)) {
                cout << fName.Data() << endl;
            }
        }
    }
}



