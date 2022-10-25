#ifndef _CONFIGMANAGER_H_
#define _CONFIGMANAGER_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

class ConfigManager {
public:
    ConfigManager();
    ~ConfigManager();
    int readConfig(char* fname);
private:
};

#endif
