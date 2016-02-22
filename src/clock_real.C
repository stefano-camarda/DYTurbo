#include "clock_real.h"

#include <sys/time.h>
#include <ctime>

double clock_real(){
    struct timeval now;
    gettimeofday(&now, NULL);
    return now.tv_sec+(now.tv_usec/1000000.0); // in sec with micro second precission
}
