#ifndef MERGER
#define MERGER

/**
 * @file merger.C
 * Description of this file
 *
 * @brief A brief description
 *
 * @author cuto <Jakub.Cuth@cern.ch>
 * @author Stefano <Stefano.Camarda@cern.ch>
 * @date 2015-11-18
 */

class OutlierRemoval{
    public :
        Merger(){};
        ~Merger(){};

        // Public methods
    private :
        // Private methods

        // Data memebers

};


void help(const char * prog){
      cout << "usage: " << prog << " <output> <input list>" << endl;
}

/**
 * Description of main program
 *
 */
int main(int argc, const char * argv[]){

    if (argc < 4)
    {
        printf("Not enough arguments (at least 1 output and 2 inputs )\n");
        help(argv[0]);
        return 1;
    }
    OutlierRemoval merger;
    //First argument is the output file
    merger.SetOutputFile(argv[1]);
    //Next arguments are input files
    for (int i = 2; i < /*argc*/ 10; i++){
        merger.AddInputFile(argv[i]);
    }
    merger.Merge();
    merger.Write();

    return 0;
}


#endif // MERGER
