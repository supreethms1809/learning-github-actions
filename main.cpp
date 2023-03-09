#include <iostream>
#include <cstdlib>
#include <ctime>
#include "CDF97.h"


using namespace std;

int main()
{
    // Set the random seed
    srand(time(nullptr));
    sperr::dim Size;
    Size.x = 8;
    Size.y = 8;
    Size.z = 8;

    // Define the dimensions of the 3D array


    // Allocate memory for the 3D array
    float*** arr3D  = new float**[Size.x];;
    float*** arr3D_new = new float**[Size.x];;
    sperr::allocate3D<float***>(arr3D, Size);
    sperr::allocate3D<float***>(arr3D_new, Size);

    // Initialize the 3D array with random numbers
    sperr::random<float***>(arr3D, Size, false);
    sperr::random<float***>(arr3D_new, Size, true);
    std::cout << "Testing the github action workflow "<<std::endl;

    sperr::print_array<float***>(arr3D,Size);

    // Free the memory allocated for the 3D array
    sperr::deallocate<float***>(arr3D,Size.x,Size.y);
    sperr::deallocate<float***>(arr3D_new,Size.x,Size.y);
    delete[] arr3D;
    delete[] arr3D_new;
    return 0;
}
