/* UserDefinedSettings' manages the parameters that have been set by the users or have predefined default values. */


// singleton pattern

#include <string>

#ifndef USER_DEFINED_SETTINGS_H
#define USER_DEFINED_SETTINGS_H

class UserDefinedSettings {

public:
   static UserDefinedSettings* Instance();
	std::string referenceFilename;
	std::string pindelConfigFilename;
	std::string pindelFilename;
	std::string bamConfigFilename;
private:
   UserDefinedSettings();
	static UserDefinedSettings* m_instance;
};

#endif
