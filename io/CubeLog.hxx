#include <iostream>
#include <iomanip>

#ifndef CUBE_LOG_LEVEL
#define CUBE_LOG_LEVEL 1
#endif

#ifndef CUBE_LOG
#define CUBE_LOG(level) if ((level) <= (CUBE_LOG_LEVEL)) std::cout << std::setfill('%') << std::setw(level+2) << " " << std::setfill(' ')
#endif

#ifndef CUBE_DEBUG_LEVEL
#define CUBE_DEBUG_LEVEL 1
#endif

#ifndef CUBE_DEBUG
#define CUBE_DEBUG(level) if ((level) <= (CUBE_DEBUG_LEVEL)) std::cout << std::setfill('#') << std::setw(level+2) << " " << std::setfill(' ') << "(" << __LINE__ << "):"
#endif

#ifndef CUBE_ERROR
#define CUBE_ERROR (std::cout <<__FILE__<<":: " << __LINE__ << ": " )
#endif
