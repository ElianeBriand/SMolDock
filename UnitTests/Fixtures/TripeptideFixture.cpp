//
// Created by eliane on 15/12/2019.
//

#include "TripeptideFixture.hpp"

TripeptideFixture::TripeptideFixture() {

    this-> tripeptide_PDB_block = "ATOM      1  N   LEU A 353      10.173  29.114  20.995  1.00 26.37           N\n"
               "ATOM      2  CA  LEU A 353      11.212  28.221  20.486  1.00 27.75           C\n"
               "ATOM      3  C   LEU A 353      10.709  26.822  20.271  1.00 30.31           C\n"
               "ATOM      4  O   LEU A 353      11.133  26.165  19.316  1.00 30.45           O\n"
               "ATOM      5  CB  LEU A 353      12.418  28.143  21.455  1.00 27.35           C\n"
               "ATOM      6  CG  LEU A 353      13.515  29.172  21.241  1.00 30.47           C\n"
               "ATOM      7  CD1 LEU A 353      14.611  28.989  22.313  1.00 24.80           C\n"
               "ATOM      8  CD2 LEU A 353      14.086  29.104  19.772  1.00 28.17           C\n"
               "ATOM      9  N   SER A 354       9.847  26.352  21.182  1.00 23.26           N\n"
               "ATOM     10  CA  SER A 354       9.436  24.954  21.193  1.00 27.45           C\n"
               "ATOM     11  C   SER A 354       8.706  24.592  19.900  1.00 30.13           C\n"
               "ATOM     12  O   SER A 354       8.867  23.486  19.379  1.00 27.28           O\n"
               "ATOM     13  CB  SER A 354       8.521  24.635  22.372  1.00 25.47           C\n"
               "ATOM     14  OG  SER A 354       7.269  25.291  22.254  1.00 26.67           O\n"
               "ATOM     15  N   GLY A 355       7.897  25.535  19.426  1.00 26.54           N\n"
               "ATOM     16  CA  GLY A 355       6.979  25.285  18.321  1.00 31.56           C\n"
               "ATOM     17  C   GLY A 355       5.788  24.428  18.724  1.00 32.14           C\n"
               "ATOM     18  O   GLY A 355       5.030  23.991  17.866  1.00 34.93           O\n"
               "ATOM     19  N   TYR A 356       5.602  24.194  20.019  1.00 30.09           N\n"
               "ATOM     20  CA  TYR A 356       4.545  23.298  20.495  1.00 31.97           C\n"
               "ATOM     21  C   TYR A 356       3.174  23.941  20.501  1.00 33.46           C\n"
               "ATOM     22  O   TYR A 356       3.032  25.134  20.775  1.00 33.20           O\n"
               "ATOM     23  CB  TYR A 356       4.813  22.795  21.926  1.00 33.21           C\n"
               "ATOM     24  CG  TYR A 356       6.115  22.070  22.134  1.00 31.28           C\n"
               "ATOM     25  CD1 TYR A 356       6.868  21.615  21.053  1.00 30.99           C\n"
               "ATOM     26  CD2 TYR A 356       6.608  21.859  23.424  1.00 27.88           C\n"
               "ATOM     27  CE1 TYR A 356       8.067  20.953  21.253  1.00 31.50           C\n"
               "ATOM     28  CE2 TYR A 356       7.800  21.208  23.627  1.00 32.71           C\n"
               "ATOM     29  CZ  TYR A 356       8.526  20.767  22.553  1.00 29.59           C\n"
               "ATOM     30  OH  TYR A 356       9.720  20.132  22.747  1.00 32.38           O\n";

    this->tripeptide_Protein = std::make_shared<sd::Protein>();
    this->tripeptide_Protein->populateFromPDBString(this-> tripeptide_PDB_block);

    this->tripeptide_iProtein = std::make_shared<sd::iProtein>();
    *(this->tripeptide_iProtein) = this->tripeptide_Protein->getiProtein();

}